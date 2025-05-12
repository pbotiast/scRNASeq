# -*- coding: utf-8 -*-
"""
Streamlit app para análisis interactivo de datos de Single-Cell RNA-seq con múltiples muestras.
Permite cargar datos 10x, realizar QC, normalización, PCA, UMAP, clustering, DEA y explorar genes.
Versión mejorada con optimización de memoria, manejo de errores y UX mejorada,
incluyendo UMAP 3D, personalización de plots y guardado/carga de parámetros.
"""
import streamlit as st
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import tempfile
import io
import traceback
import plotly.express as px
import plotly.graph_objects as go 
import difflib
import gzip
from scipy import sparse
from scipy.stats import zscore 
import json 
import markdown

# Configuración de la página
st.set_page_config(layout="wide", initial_sidebar_state="expanded")
st.title("Analizador Interactivo de Single-Cell RNA-seq Avanzado")

# --- Funciones Auxiliares ---
def validate_10x_files(matrix_file, features_file, barcodes_file):
    files_to_check = {
        "Matrix": (matrix_file, "%%MatrixMarket"),
        "Features": (features_file, None), 
        "Barcodes": (barcodes_file, None)
    }
    try:
        for name, (file_obj, expected_start) in files_to_check.items():
            if file_obj is None: 
                raise ValueError(f"Archivo {name} no proporcionado.")
            
            content = file_obj.getbuffer()
            is_gz = file_obj.name.endswith(".gz")
            
            if is_gz:
                try:
                    content = gzip.decompress(content)
                except gzip.BadGzipFile:
                    raise ValueError(f"Archivo {name} ({file_obj.name}) parece .gz pero no se pudo descomprimir.")
            
            if expected_start: 
                with io.BytesIO(content) as f_content:
                    first_line = f_content.readline().decode(errors='ignore')
                    if not first_line.startswith(expected_start):
                        raise ValueError(f"Archivo {name} ({file_obj.name}) no es válido. Esperaba '{expected_start}', obtuvo '{first_line[:50]}...'")
            else: 
                with io.BytesIO(content) as f_content:
                    pd.read_csv(f_content, sep="\t", header=None, nrows=5, comment='%', on_bad_lines='warn') 
        return True
    except ValueError as ve: 
        st.error(f"Error de validación: {ve}")
        return False
    except Exception as e: 
        st.error(f"Error procesando archivo para validación: {e}")
        return False

def load_10x_data(matrix_file, features_file, barcodes_file, sample_name):
    with tempfile.TemporaryDirectory() as temp_dir:
        file_map = {
            "matrix.mtx": matrix_file,
            "features.tsv": features_file, 
            "barcodes.tsv": barcodes_file
        }
        for base_name, uploaded_file in file_map.items():
            filename_in_temp = base_name
            if uploaded_file.name.endswith(".gz"):
                filename_in_temp += ".gz"
            
            with open(os.path.join(temp_dir, filename_in_temp), "wb") as f:
                f.write(uploaded_file.getbuffer()) 
        
        adata = sc.read_10x_mtx(temp_dir, var_names='gene_symbols', cache=False)
        adata.var_names_make_unique()
        adata.obs['sample'] = sample_name
        adata.X = sparse.csr_matrix(adata.X)
    return adata

def fig_to_bytes(fig, format='png'):
    buf = io.BytesIO()
    fig.savefig(buf, format=format, bbox_inches='tight', dpi=300)
    buf.seek(0)
    return buf.getvalue()

def suggest_genes(input_genes, valid_genes_lower_map):
    suggestions = {}
    valid_genes_lower_list = list(valid_genes_lower_map.keys())    

    for gene_input_raw in input_genes:
        gene_input_lower = gene_input_raw.lower()
        if gene_input_lower not in valid_genes_lower_list:
            matches_lower = difflib.get_close_matches(gene_input_lower, valid_genes_lower_list, n=3, cutoff=0.6)
            if matches_lower:
                suggestions[gene_input_raw] = [valid_genes_lower_map[m] for m in matches_lower]
    return suggestions

# --- Inicialización de st.session_state ---
default_values = {
    "num_samples": 1, "sample_files": {}, "min_genes": 200, "min_cells": 3, "mito_prefix": "MT-",
    "max_mito_pct": 15, "n_top_genes_hvg": 2000, "n_pcs": 30, "n_neighbors_val": 15, "leiden_res": 0.8,
    "n_top_markers": 10, "adata_raw": None, "adata_processed": None, "adata_hvg_subset": None,
    "analysis_done": False, "marker_genes_df": None, "condition_assignments": {},
    "dea_group1": None, "dea_group2": None, "dea_cluster_scope": "Todos los Clústeres",
    "dea_n_genes_display": 25, "dea_lfc_cutoff": 0.5, "dea_pval_cutoff": 0.05,
    "dea_results_df": None, "dea_comparison_str": "", "gene_explorer_input": "",
    "leiden_flavor": "igraph", 
    "umap_init_pos": "random", "umap_n_neighbors": 15, "umap_min_dist": 0.5, "calc_umap_3d": False,
    "plot_palette": "viridis", "plot_point_size": 20, 
    "heatmap_top_n_genes": 3 
}
for key, value in default_values.items():
    if key not in st.session_state:
        st.session_state[key] = value

# --- Sidebar ---
with st.sidebar:
    st.image("https://www.biogenouest.org/wp-content/uploads/2020/02/Logo-ScanPy-195x150.png", width=100) 
    st.header("Configuración del Análisis")
    
    with st.expander("Guardar/Cargar Configuración", expanded=False):
        if st.button("Guardar Configuración Actual", key="save_params_btn_sidebar"):
            params_to_save = {k: st.session_state[k] for k in default_values if not k.startswith("adata") and k != "sample_files"} 
            params_to_save.pop("analysis_done", None) 
            params_to_save.pop("marker_genes_df", None)
            params_to_save.pop("dea_results_df", None)
            
            params_json = json.dumps(params_to_save, indent=4)
            st.download_button(
                label="Descargar archivo de parámetros (.json)",
                data=params_json,
                file_name=f"scRNAseq_app_params_{pd.Timestamp.now().strftime('%Y%m%d')}.json",
                mime="application/json",
                key="download_params_json_sidebar_btn"
            )

        uploaded_params_file_sidebar = st.file_uploader("Cargar Configuración (.json)", type="json", key="upload_params_sidebar")
        if uploaded_params_file_sidebar is not None:
            try:
                loaded_params = json.load(uploaded_params_file_sidebar)
                for key, value_loaded in loaded_params.items():
                    if key in st.session_state:
                        st.session_state[key] = value_loaded
                st.success("Parámetros cargados. La interfaz se ha actualizado.")
                st.info("Por favor, revisa los parámetros y vuelve a ejecutar el pipeline si es necesario.")
            except Exception as e_load_params_sidebar:
                st.error(f"Error al cargar parámetros: {e_load_params_sidebar}")

    with st.expander("1. Carga de Datos", expanded=True):
        st.session_state.num_samples = st.number_input(
            "Número de muestras", min_value=1, max_value=10, value=st.session_state.num_samples, step=1, key="num_samples_main"
        )
        for i in range(st.session_state.num_samples):
            sample_widget_key_prefix = f"sample_{i}"
            st.subheader(f"Muestra {i+1}")
            
            s_name_key = f"sample_name_{i}"
            if s_name_key not in st.session_state.sample_files:
                 st.session_state.sample_files[s_name_key] = f"Muestra{i+1}"

            st.session_state.sample_files[s_name_key] = st.text_input(
                f"Nombre Muestra {i+1}", value=st.session_state.sample_files[s_name_key], key=f"{sample_widget_key_prefix}_name"
            )
            for file_type_key, label in zip(["matrix_file", "features_file", "barcodes_file"], 
                                            ["Matrix (.mtx/.mtx.gz)", "Features (.tsv/.tsv.gz)", "Barcodes (.tsv/.tsv.gz)"]):
                full_key = f"{file_type_key}_{i}"
                if full_key not in st.session_state.sample_files:
                    st.session_state.sample_files[full_key] = None
                st.session_state.sample_files[full_key] = st.file_uploader(
                    f"{label}", type=["mtx", "tsv", "gz"], key=f"{sample_widget_key_prefix}_{file_type_key}"
                )
            if i < st.session_state.num_samples - 1:
                st.markdown("---")

    with st.expander("2. Parámetros del Pipeline"):
        st.session_state.min_genes = st.slider("Mínimo genes/célula", 50, 1000, st.session_state.min_genes, key="min_genes_slider", help="Filtra células con menos de N genes detectados.")
        st.session_state.min_cells = st.slider("Mínimo células/gen", 1, 50, st.session_state.min_cells, key="min_cells_slider", help="Filtra genes detectados en menos de N células.")
        st.session_state.mito_prefix = st.text_input("Prefijo genes mitocondriales", st.session_state.mito_prefix, key="mito_prefix_input", help="Ej: 'MT-' para humano, 'mt-' para ratón.")
        st.session_state.max_mito_pct = st.slider("Máx % cuentas mitocondriales", 1, 100, st.session_state.max_mito_pct, key="max_mito_slider", help="Filtra células con alto % de cuentas mitocondriales.")
        st.session_state.n_top_genes_hvg = st.slider("Nº HVGs a seleccionar", 500, 5000, st.session_state.n_top_genes_hvg, key="n_hvg_slider")
        st.session_state.n_pcs = st.slider("Nº PCs (para PCA y Vecinos)", 5, 100, st.session_state.n_pcs, key="n_pcs_slider", help="Mínimo 5 recomendado para 'arpack'.")
        st.session_state.n_neighbors_val = st.slider("Nº Vecinos (para grafo KNN)", 2, 50, st.session_state.n_neighbors_val, key="n_neighbors_slider", help="Usado para construir el grafo de vecinos para UMAP y Leiden.")
        st.session_state.leiden_res = st.slider("Resolución Leiden", 0.1, 2.5, st.session_state.leiden_res, 0.1, key="leiden_res_slider", help="Mayor resolución = más clústeres.")
        st.session_state.n_top_markers = st.slider("Nº marcadores a mostrar/clúster", 1, 25, st.session_state.n_top_markers, key="n_markers_slider")
        
        leiden_flavors_list = ["igraph", "leidenalg"]
        st.session_state.leiden_flavor = st.selectbox("Backend Leiden", leiden_flavors_list, 
                                                      index=leiden_flavors_list.index(st.session_state.leiden_flavor) if st.session_state.leiden_flavor in leiden_flavors_list else 0, 
                                                      key="leiden_flavor_select", help="'igraph' es generalmente recomendado.")
        
        st.subheader("Parámetros UMAP")
        st.session_state.calc_umap_3d = st.checkbox("Calcular también UMAP 3D", value=st.session_state.calc_umap_3d, key="calc_3d_umap_check")
        umap_init_options_list = ["spectral", "random", "pca"]
        st.session_state.umap_init_pos = st.selectbox("Inicialización UMAP", umap_init_options_list, 
                                                      index=umap_init_options_list.index(st.session_state.umap_init_pos) if st.session_state.umap_init_pos in umap_init_options_list else 1, 
                                                      key="umap_init_select", help="'random' puede ser más estable con algunas versiones de NumPy.")
        st.session_state.umap_n_neighbors = st.slider("Nº Vecinos UMAP (para embedding)", 2, 200, st.session_state.umap_n_neighbors, key="umap_n_neighbors_slider", help="Controla el balance entre estructura local y global en UMAP.")
        st.session_state.umap_min_dist = st.slider("Distancia Mínima UMAP", 0.0, 1.0, st.session_state.umap_min_dist, 0.01, key="umap_min_dist_slider", help="Controla cuán agrupados estarán los puntos en UMAP.")

    with st.expander("Personalización de Plots", expanded=False):
        palettes_options = ["tab10", "tab20", "Set3", "Paired", "viridis", "plasma", "magma", "cividis"]
        st.session_state.plot_palette = st.selectbox("Paleta de Colores (Clusters/Muestras)", palettes_options, 
                                                     index=palettes_options.index(st.session_state.plot_palette) if st.session_state.plot_palette in palettes_options else 0,
                                                     key="plot_palette_select")
        st.session_state.plot_point_size = st.slider("Tamaño de Puntos UMAP (aprox.)", 10, 150, st.session_state.plot_point_size, 5, key="plot_point_size_slider", help="Valor de 's' para sc.pl.umap (multiplicado por un factor).")
        st.session_state.heatmap_top_n_genes = st.slider("Nº genes/clúster para Heatmap Marcadores", 1, 10, st.session_state.heatmap_top_n_genes, key="heatmap_genes_slider")

    if st.session_state.analysis_done and st.session_state.adata_processed is not None:
        with st.expander("3. Análisis Diferencial (DEA)"):
            adata_for_dea_config = st.session_state.adata_processed
            samples_in_adata = sorted(adata_for_dea_config.obs['sample'].unique().tolist())
            st.subheader("Asignar Condiciones")
            
            current_assignments_dea = {s: st.session_state.condition_assignments.get(s, f"Cond_{s.replace(' ','_')}") for s in samples_in_adata}
            for sample_name_dea_config in samples_in_adata:
                current_assignments_dea[sample_name_dea_config] = st.text_input(
                    f"Condición para {sample_name_dea_config}", 
                    value=current_assignments_dea[sample_name_dea_config], 
                    key=f"cond_assign_{sample_name_dea_config}"
                )
            st.session_state.condition_assignments = current_assignments_dea

            unique_defined_conditions = sorted(list(set(c for c in st.session_state.condition_assignments.values() if c and c.strip())))

            if len(unique_defined_conditions) >= 2:
                st.subheader("Seleccionar Grupos para Comparación")
                col1_dea, col2_dea = st.columns(2)
                with col1_dea:
                    st.session_state.dea_group1 = st.selectbox("Grupo 1 (Referencia)", unique_defined_conditions, 
                                                               index=unique_defined_conditions.index(st.session_state.dea_group1) if st.session_state.dea_group1 in unique_defined_conditions else 0,
                                                               key="dea_g1_select")
                with col2_dea:
                    options_g2 = [c for c in unique_defined_conditions if c != st.session_state.dea_group1]
                    if not options_g2 and unique_defined_conditions:
                        options_g2 = [c for c in unique_defined_conditions if c != st.session_state.dea_group1] or (unique_defined_conditions[1:] if len(unique_defined_conditions)>1 else unique_defined_conditions)
                    
                    st.session_state.dea_group2 = st.selectbox("Grupo 2 (Comparación)", options_g2, 
                                                               index=options_g2.index(st.session_state.dea_group2) if st.session_state.dea_group2 in options_g2 and options_g2 else 0,
                                                               key="dea_g2_select")

                clusters_for_dea = ["Todos los Clústeres"] + sorted(adata_for_dea_config.obs['leiden_clusters'].astype(str).unique().tolist())
                st.session_state.dea_cluster_scope = st.selectbox("Ámbito DEA", clusters_for_dea, 
                                                                  index=clusters_for_dea.index(st.session_state.dea_cluster_scope) if st.session_state.dea_cluster_scope in clusters_for_dea else 0,
                                                                  key="dea_scope_select")
                
                st.session_state.dea_n_genes_display = st.slider("Nº genes DEA a mostrar", 10, 200, st.session_state.dea_n_genes_display, key="dea_ngenes_slider")
                st.session_state.dea_lfc_cutoff = st.number_input("Log2FC cutoff (para Volcano)", 0.0, value=st.session_state.dea_lfc_cutoff, step=0.1, key="dea_lfc_input")
                st.session_state.dea_pval_cutoff = st.number_input("P-adj cutoff (para Volcano)", 0.0, 1.0, value=st.session_state.dea_pval_cutoff, step=0.01, format="%.3f", key="dea_pval_input")

                adata_for_dea_preview = adata_for_dea_config.copy()
                adata_for_dea_preview.obs['condition_temp_dea'] = adata_for_dea_preview.obs['sample'].map(st.session_state.condition_assignments)
                
                valid_cells_for_preview = adata_for_dea_preview.obs['condition_temp_dea'].notna()
                if valid_cells_for_preview.any():
                    counts_df = adata_for_dea_preview[valid_cells_for_preview].obs.groupby(['condition_temp_dea', 'leiden_clusters'], observed=True).size().unstack(fill_value=0)
                    st.write("**Conteos de Células por Condición y Clúster (para condiciones asignadas):**")
                    st.dataframe(counts_df)
                else:
                    st.warning("No hay células con condiciones asignadas para mostrar conteos.")

                if st.button("Ejecutar DEA", key="run_dea_btn"):
                    if not st.session_state.dea_group1 or not st.session_state.dea_group2:
                        st.error("Por favor, selecciona ambos grupos para la comparación.")
                    elif st.session_state.dea_group1 == st.session_state.dea_group2:
                        st.error("Los grupos de comparación deben ser diferentes.")
                    else:
                        with st.spinner("Ejecutando DEA..."):
                            try:
                                adata_filtered_for_dea = adata_for_dea_preview[
                                    adata_for_dea_preview.obs['condition_temp_dea'].isin([st.session_state.dea_group1, st.session_state.dea_group2]) &
                                    adata_for_dea_preview.obs['condition_temp_dea'].notna() 
                                ].copy()
                                scope_msg_dea = ""
                                if st.session_state.dea_cluster_scope != "Todos los Clústeres":
                                    scope_msg_dea = f" (Clúster {st.session_state.dea_cluster_scope})"
                                    adata_filtered_for_dea = adata_filtered_for_dea[
                                        adata_filtered_for_dea.obs['leiden_clusters'] == st.session_state.dea_cluster_scope
                                    ].copy()
                                cell_counts_in_groups = adata_filtered_for_dea.obs['condition_temp_dea'].value_counts()
                                if (st.session_state.dea_group1 not in cell_counts_in_groups or cell_counts_in_groups[st.session_state.dea_group1] < 3 or
                                    st.session_state.dea_group2 not in cell_counts_in_groups or cell_counts_in_groups[st.session_state.dea_group2] < 3):
                                    st.error(f"Insuficientes células para DEA en los grupos seleccionados{scope_msg_dea}. Conteos: {cell_counts_in_groups.to_dict()}")
                                    st.session_state.dea_results_df = None
                                else:
                                    sc.tl.rank_genes_groups(
                                        adata_filtered_for_dea, 'condition_temp_dea', groups=[st.session_state.dea_group2],
                                        reference=st.session_state.dea_group1, method='wilcoxon', key_added='rank_genes_dea', use_raw=False
                                    )
                                    res_dea_uns = adata_filtered_for_dea.uns['rank_genes_dea']
                                    group_name_in_uns = st.session_state.dea_group2
                                    df_dea_results = pd.DataFrame({
                                        'Gene': res_dea_uns['names'][group_name_in_uns],
                                        'Log2FC': res_dea_uns['logfoldchanges'][group_name_in_uns],
                                        'P-Value': res_dea_uns['pvals'][group_name_in_uns],
                                        'P-Value Adj': res_dea_uns['pvals_adj'][group_name_in_uns],
                                        'Scores': res_dea_uns['scores'][group_name_in_uns]
                                    })
                                    st.session_state.dea_results_df = df_dea_results.sort_values('P-Value Adj')
                                    st.session_state.dea_comparison_str = f"{st.session_state.dea_group2} vs {st.session_state.dea_group1}{scope_msg_dea}"
                                    st.success(f"DEA completado para {st.session_state.dea_comparison_str}.")
                            except Exception as e_dea_run:
                                st.error(f"Error en DEA: {e_dea_run}\n{traceback.format_exc()}")
                                st.session_state.dea_results_df = None
            else:
                st.info("Define al menos dos condiciones válidas para habilitar el DEA.")

    all_files_provided = True
    if st.session_state.num_samples > 0:
        for i in range(st.session_state.num_samples):
            if not (st.session_state.sample_files.get(f"matrix_file_{i}") and \
                    st.session_state.sample_files.get(f"features_file_{i}") and \
                    st.session_state.sample_files.get(f"barcodes_file_{i}")):
                all_files_provided = False
                break
    else:
        all_files_provided = False

    if all_files_provided:
        if st.button("Cargar y Concatenar Datos", key="load_concat_btn_main", type="primary"):
            validation_messages = []
            all_valid_globally = True
            for i in range(st.session_state.num_samples):
                sample_name_val = st.session_state.sample_files.get(f"sample_name_{i}", f"Muestra {i+1}")
                matrix_f_val = st.session_state.sample_files.get(f"matrix_file_{i}")
                features_f_val = st.session_state.sample_files.get(f"features_file_{i}")
                barcodes_f_val = st.session_state.sample_files.get(f"barcodes_file_{i}")

                if not (matrix_f_val and features_f_val and barcodes_f_val):
                    all_valid_globally = False
                    validation_messages.append(f"Muestra {i+1} ({sample_name_val}): Faltan archivos.")
                    continue 

                validation_messages.append(f"Validando Muestra {i+1} ({sample_name_val}):")
                if not validate_10x_files(matrix_f_val, features_f_val, barcodes_f_val):
                    all_valid_globally = False
                    validation_messages.append(f" -> Archivos de Muestra {i+1} ({sample_name_val}) NO son válidos.")
                else:
                    validation_messages.append(f" -> Archivos de Muestra {i+1} ({sample_name_val}) VALIDADOS.")
            
            with st.expander("Registro de Validación de Archivos", expanded=not all_valid_globally):
                for msg in validation_messages:
                    if "NO son válidos" in msg or "Faltan archivos" in msg : st.error(msg)
                    else: st.info(msg)

            if all_valid_globally:
                with st.spinner("Cargando y concatenando datos..."):
                    try:
                        adatas_dict = {} 
                        for i in range(st.session_state.num_samples):
                            user_sample_name = st.session_state.sample_files[f"sample_name_{i}"]
                            adata_sample = load_10x_data(
                                st.session_state.sample_files[f"matrix_file_{i}"],
                                st.session_state.sample_files[f"features_file_{i}"],
                                st.session_state.sample_files[f"barcodes_file_{i}"],
                                user_sample_name 
                            )
                            if adata_sample is None:
                                raise ValueError(f"Fallo al cargar Muestra {i+1} ({user_sample_name}).")
                            adatas_dict[user_sample_name] = adata_sample
                        
                        st.session_state.adata_raw = ad.concat(adatas_dict, label='sample', index_unique='-', join='outer', fill_value=0)
                        
                        st.session_state.adata_processed = None 
                        st.session_state.analysis_done = False
                        for key_to_reset in ["marker_genes_df", "dea_results_df", "dea_comparison_str", "dea_group1", "dea_group2", "dea_cluster_scope"]:
                            st.session_state[key_to_reset] = default_values.get(key_to_reset)

                        st.success(f"Carga completada: {st.session_state.adata_raw.n_obs} células, {st.session_state.adata_raw.n_vars} genes.")
                        print("DEBUG: Muestras en adata_raw.obs['sample']:", st.session_state.adata_raw.obs['sample'].unique())
                    except Exception as e_load:
                        st.error(f"Error durante la carga: {e_load}")
                        st.error(traceback.format_exc())
                        st.session_state.adata_raw = None

    elif st.session_state.num_samples > 0 :
        st.warning("Por favor, sube todos los archivos para cada muestra para habilitar la carga.")

    if st.session_state.adata_raw is not None:
        if st.button("Ejecutar Pipeline Principal", key="run_pipeline_btn_main_v2", type="primary"):
            adata = st.session_state.adata_raw.copy() 
            st.session_state.adata_hvg_subset = None 
            
            with st.spinner("Ejecutando pipeline... Esto puede tardar varios minutos."):
                progress_bar = st.progress(0)
                status_text = st.empty()
                try:
                    # 1. QC
                    status_text.text("Paso 1/8: Control de Calidad (QC)...")
                    sc.pp.filter_cells(adata, min_genes=st.session_state.min_genes)
                    sc.pp.filter_genes(adata, min_cells=st.session_state.min_cells)
                    if adata.n_obs == 0 or adata.n_vars == 0: raise ValueError(f"QC inicial eliminó todas las células ({adata.n_obs}) o genes ({adata.n_vars}).")
                    
                    adata.var['mt'] = adata.var_names.str.startswith(st.session_state.mito_prefix)
                    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
                    adata = adata[adata.obs.pct_counts_mt < st.session_state.max_mito_pct, :].copy()
                    if adata.n_obs == 0: raise ValueError("QC de % mitocondrial eliminó todas las células.")
                    progress_bar.progress(12)

                    # 2. HVGs
                    status_text.text("Paso 2/8: Selección de Genes Altamente Variables (HVGs)...")
                    print(f"DEBUG HVG: Shape antes de HVG: {adata.shape}, batch_key='sample'")
                    sc.pp.highly_variable_genes(adata, n_top_genes=st.session_state.n_top_genes_hvg, flavor='seurat_v3', batch_key='sample')
                    if 'highly_variable' not in adata.var.columns or adata.var['highly_variable'].sum() < 10:
                        raise ValueError(f"Muy pocos HVGs detectados ({adata.var.get('highly_variable', pd.Series(False)).sum()}).")
                    progress_bar.progress(25)

                    # 3. Normalización y Log-transformación
                    status_text.text("Paso 3/8: Normalización y Transformación Logarítmica...")
                    sc.pp.normalize_total(adata, target_sum=1e4)
                    sc.pp.log1p(adata)
                    progress_bar.progress(37)

                    # 4. Crear subconjunto para HVG
                    status_text.text("Paso 4/8: Creando subconjunto de HVGs...")
                    st.session_state.adata_hvg_subset = adata[:, adata.var.highly_variable].copy()
                    if st.session_state.adata_hvg_subset.n_vars == 0: raise ValueError("Subconjunto de HVG sin variables.")
                    progress_bar.progress(50)
                    
                    # 5. Escalado
                    status_text.text("Paso 5/8: Escalado de HVGs...")
                    sc.pp.scale(st.session_state.adata_hvg_subset, max_value=10)
                    progress_bar.progress(62)

                    # 6. PCA, Vecinos
                    status_text.text("Paso 6/8: PCA y Cálculo de Vecinos...")
                    n_obs_hvg, n_vars_hvg = st.session_state.adata_hvg_subset.shape
                    print(f"DEBUG PCA: Shape adata_hvg_subset={n_obs_hvg}x{n_vars_hvg}, n_pcs_slider={st.session_state.n_pcs}")

                    current_n_pcs_val = int(st.session_state.n_pcs)
                    limit_pcs = min(n_obs_hvg, n_vars_hvg)
                    if limit_pcs <= 1: raise ValueError(f"Dimensiones para PCA ({n_obs_hvg}x{n_vars_hvg}) muy pequeñas (min_dim <=1).")
                    valid_max_n_pcs = limit_pcs - 1 
                    if current_n_pcs_val > valid_max_n_pcs : current_n_pcs_val = valid_max_n_pcs; st.warning(f"Nº PCs ajustado a {current_n_pcs_val}")
                    if current_n_pcs_val <= 0 : current_n_pcs_val = 1; st.warning(f"Nº PCs <=0, ajustado a {current_n_pcs_val}")
                    
                    print(f"DEBUG PCA: Usando n_comps={current_n_pcs_val}")
                    sc.tl.pca(st.session_state.adata_hvg_subset, svd_solver='arpack', n_comps=current_n_pcs_val, random_state=0)
                    
                    if 'X_pca' in st.session_state.adata_hvg_subset.obsm:
                        X_pca_check = st.session_state.adata_hvg_subset.obsm['X_pca']
                        if np.any(np.isnan(X_pca_check)) or np.any(np.isinf(X_pca_check)): raise ValueError("NaNs o Infs en X_pca.")
                        print("DEBUG PCA: No NaNs/Infs en X_pca.")
                    else: raise KeyError("'X_pca' no encontrado post-PCA.")

                    current_n_neighbors_val = int(st.session_state.n_neighbors_val)
                    if n_obs_hvg <= 1: raise ValueError(f"Pocas células ({n_obs_hvg}) para calcular vecinos.")
                    valid_max_n_neighbors = n_obs_hvg - 1 if n_obs_hvg >1 else 1 
                    if current_n_neighbors_val > valid_max_n_neighbors: current_n_neighbors_val = valid_max_n_neighbors; st.warning(f"Nº Vecinos ajustado a {current_n_neighbors_val}")
                    if current_n_neighbors_val <= 0: current_n_neighbors_val = 1; st.warning(f"Nº Vecinos <=0, ajustado a {current_n_neighbors_val}")

                    print(f"DEBUG Neighbors: Usando n_neighbors={current_n_neighbors_val}, n_pcs={current_n_pcs_val}")
                    sc.pp.neighbors(st.session_state.adata_hvg_subset, n_neighbors=current_n_neighbors_val, n_pcs=current_n_pcs_val, random_state=0)
                    progress_bar.progress(75)

                    # 7. UMAP
                    status_text.text("Paso 7/8: Cálculo de UMAP...")
                    print(f"DEBUG UMAP: init={st.session_state.umap_init_pos}, n_neighbors={st.session_state.umap_n_neighbors}, min_dist={st.session_state.umap_min_dist}")
                    X_pca_for_umap = st.session_state.adata_hvg_subset.obsm['X_pca']
                    
                    # UMAP 2D
                    try:
                        import umap 
                        reducer_2d = umap.UMAP(
                            n_neighbors=int(st.session_state.umap_n_neighbors), n_components=2, metric='euclidean',
                            min_dist=float(st.session_state.umap_min_dist), init=st.session_state.umap_init_pos,
                            random_state=42, verbose=False 
                        )
                        embedding_2d = reducer_2d.fit_transform(X_pca_for_umap)
                        st.session_state.adata_hvg_subset.obsm['X_umap'] = embedding_2d
                        print("DEBUG UMAP 2D: API directa completada.")
                    except Exception as e_umap_direct_api:
                        st.warning(f"UMAP 2D con API directa falló ({e_umap_direct_api}). Intentando con sc.tl.umap...")
                        try:
                            umap_init_fallback = st.session_state.umap_init_pos if st.session_state.umap_init_pos != 'pca' else 'spectral'
                            sc.tl.umap(st.session_state.adata_hvg_subset, 
                                       n_neighbors=int(st.session_state.umap_n_neighbors), 
                                       min_dist=float(st.session_state.umap_min_dist),
                                       init_pos=umap_init_fallback, 
                                       random_state=0) 
                            print("DEBUG UMAP 2D: sc.tl.umap (fallback) completado.")
                        except Exception as e_sc_umap_fallback:
                             st.error(f"Fallback sc.tl.umap para UMAP 2D también falló: {e_sc_umap_fallback}.")
                             st.session_state.adata_hvg_subset.obsm.pop('X_umap', None) # Eliminar si falló

                    # UMAP 3D (si está seleccionado)
                    if st.session_state.calc_umap_3d:
                        status_text.text("Calculando UMAP 3D...")
                        try:
                            import umap
                            reducer_3d = umap.UMAP(
                                n_neighbors=int(st.session_state.umap_n_neighbors), n_components=3, metric='euclidean',
                                min_dist=float(st.session_state.umap_min_dist), init=st.session_state.umap_init_pos,
                                random_state=42, verbose=False
                            )
                            embedding_3d = reducer_3d.fit_transform(X_pca_for_umap)
                            st.session_state.adata_hvg_subset.obsm['X_umap_3d'] = embedding_3d
                            print("DEBUG UMAP 3D: API directa completada.")
                        except Exception as e_umap3d_direct:
                            st.warning(f"UMAP 3D con API directa falló ({e_umap3d_direct}). No se generará UMAP 3D.")
                            st.session_state.adata_hvg_subset.obsm.pop('X_umap_3d', None)
                    progress_bar.progress(87)


                    # 8. Clustering y Transferencia
                    status_text.text("Paso 8/8: Clustering Leiden y finalización...")
                    print(f"DEBUG Leiden: flavor={st.session_state.leiden_flavor}, res={st.session_state.leiden_res}")
                    sc.tl.leiden(
                        st.session_state.adata_hvg_subset, resolution=st.session_state.leiden_res,
                        key_added="leiden_clusters", flavor=st.session_state.leiden_flavor,
                        n_iterations=2 if st.session_state.leiden_flavor == "igraph" else -1,
                        random_state=0, directed=False
                    )
                    adata.obs['leiden_clusters'] = st.session_state.adata_hvg_subset.obs['leiden_clusters'].astype('category')
                    
                    if 'X_umap' in st.session_state.adata_hvg_subset.obsm:
                        adata.obsm['X_umap'] = st.session_state.adata_hvg_subset.obsm['X_umap']
                    if 'X_umap_3d' in st.session_state.adata_hvg_subset.obsm: # Transferir UMAP 3D si existe
                        adata.obsm['X_umap_3d'] = st.session_state.adata_hvg_subset.obsm['X_umap_3d']
                    if 'X_pca' in st.session_state.adata_hvg_subset.obsm: 
                         adata.obsm['X_pca_hvg'] = st.session_state.adata_hvg_subset.obsm['X_pca'] 

                    # Marcadores
                    adata.obs['leiden_clusters'] = adata.obs['leiden_clusters'].astype('category')
                    sc.tl.rank_genes_groups(adata, 'leiden_clusters', method='wilcoxon', key_added='rank_genes_leiden_clusters', use_raw=False)
                    
                    res_markers = adata.uns.get('rank_genes_leiden_clusters', {})
                    marker_data_list = []
                    if res_markers and 'names' in res_markers and hasattr(res_markers['names'], 'dtype') and res_markers['names'].dtype.names is not None:
                        cluster_ids_with_names_markers = res_markers['names'].dtype.names
                        for grp_marker in cluster_ids_with_names_markers:
                            # Comprobación más robusta de la existencia de claves para este grupo en todos los arrays
                            if all(grp_marker in res_markers.get(field, {}).dtype.names for field in ['names', 'scores', 'logfoldchanges', 'pvals_adj'] if hasattr(res_markers.get(field), 'dtype')):
                                num_avail_markers = len(res_markers['names'][grp_marker])
                                markers_to_fetch = min(st.session_state.n_top_markers, num_avail_markers)
                                for i_marker in range(markers_to_fetch):
                                    # Comprobación de longitud individual (aunque deberían ser iguales)
                                    if all(i_marker < len(res_markers[field][grp_marker]) for field in ['names', 'scores', 'logfoldchanges', 'pvals_adj']):
                                        marker_data_list.append({
                                            'Cluster': grp_marker, 'Rank': i_marker + 1, 
                                            'Gene': res_markers['names'][grp_marker][i_marker], 
                                            'Score': res_markers['scores'][grp_marker][i_marker],
                                            'Log2FC': res_markers['logfoldchanges'][grp_marker][i_marker], 
                                            'P-Value Adj': res_markers['pvals_adj'][grp_marker][i_marker]
                                        })
                                    else: print(f"Warn: Discrepancia longitud arrays marcadores para clúster '{grp_marker}', índice {i_marker}")
                            else: print(f"Warn: Faltan campos para marcadores del clúster '{grp_marker}'")
                    st.session_state.marker_genes_df = pd.DataFrame(marker_data_list) if marker_data_list else pd.DataFrame()
                    
                    st.session_state.adata_processed = adata 
                    st.session_state.analysis_done = True
                    progress_bar.progress(100)
                    status_text.empty()
                    st.balloons()

                except ValueError as ve_pipeline:
                    st.error(f"Error de Valor en Pipeline: {ve_pipeline}")
                    st.error(traceback.format_exc())
                except KeyError as ke_pipeline:
                    st.error(f"Error de Clave en Pipeline: {ke_pipeline}")
                    st.error(traceback.format_exc())
                except Exception as e_pipeline:
                    st.error(f"Error Inesperado en Pipeline: {e_pipeline}")
                    st.error(traceback.format_exc())
                finally: 
                    if 'analysis_done' not in st.session_state or not st.session_state.analysis_done:
                        st.session_state.analysis_done = False
                        st.session_state.adata_processed = None
                        st.session_state.adata_hvg_subset = None
                        for key_to_reset in default_values: 
                            if not key_to_reset.startswith("adata") and key_to_reset not in ["sample_files", "num_samples", "analysis_done"]:
                                st.session_state[key_to_reset] = default_values[key_to_reset]

# --- Sección de Resultados ---
st.markdown("---")
st.header("Resultados del Análisis")
st.subheader("🔬 Explorador de Expresión Génica")
st.session_state.gene_explorer_input = st.text_area(
    "Ingresa nombres de genes (separados por coma, espacio o nueva línea):", 
    value=st.session_state.gene_explorer_input, 
    key="gene_explorer_main_input_v3", 
    height=100
)

if st.session_state.analysis_done and st.session_state.adata_processed is not None:
    adata_display = st.session_state.adata_processed 
    
    valid_genes_lower_map_display = {g.lower(): g for g in adata_display.var_names}

    try:
        with tempfile.NamedTemporaryFile(delete=False, suffix=".h5ad") as tmp_h5ad_file:
            filepath_h5ad = tmp_h5ad_file.name
        adata_display.write_h5ad(filepath_h5ad)

        with open(filepath_h5ad, "rb") as f_h5ad_read:
            st.sidebar.download_button(
                "Descargar AnnData Procesado (.h5ad)", 
                f_h5ad_read.read(), 
                f"processed_adata_scRNAseq_{pd.Timestamp.now().strftime('%Y%m%d_%H%M')}.h5ad",
                "application/octet-stream",
                key="download_adata_button_sidebar_v2" 
            )
        os.remove(filepath_h5ad) 
    except Exception as e_dl_sidebar: 
        st.sidebar.error(f"Error descarga AnnData: {e_dl_sidebar}")
        if 'filepath_h5ad' in locals() and os.path.exists(filepath_h5ad):
            try: os.remove(filepath_h5ad)
            except: pass

    genes_input_list_raw = [g.strip() for g in st.session_state.gene_explorer_input.replace(',', ' ').replace('\n', ' ').split() if g.strip()]
    genes_to_visualize_list = []
    genes_not_found_in_adata_list = []
    for gene_raw in genes_input_list_raw:
        gene_lower = gene_raw.lower()
        if gene_lower in valid_genes_lower_map_display:
            genes_to_visualize_list.append(valid_genes_lower_map_display[gene_lower])
        else:
            genes_not_found_in_adata_list.append(gene_raw)
    genes_to_visualize_list = list(dict.fromkeys(genes_to_visualize_list))

    if genes_not_found_in_adata_list:
        gene_suggestions = suggest_genes(genes_not_found_in_adata_list, valid_genes_lower_map_display)
        warning_msg = f"Genes no encontrados: {', '.join(genes_not_found_in_adata_list)}."
        if gene_suggestions:
            suggestions_str_list = [f"{gene_orig} -> (quizás: {', '.join(sugg_list)}?)" for gene_orig, sugg_list in gene_suggestions.items()]
            warning_msg += f" Sugerencias: {'; '.join(suggestions_str_list)}"
        st.warning(warning_msg)

    tab_titles_main = ["📊 UMAPs", "🔬 Marcadores", "🔥 Heatmap Marcadores", "🧬 QC", "📈 DEA", "🧬 Explorador Genes", "ℹ️ Info"]
    tab_umaps, tab_markers, tab_heatmap, tab_qc, tab_dea, tab_gene_explorer, tab_info = st.tabs(tab_titles_main)
    
    sc.set_figure_params(vector_friendly=True, format='png', dpi_save=300, color_map=st.session_state.plot_palette)

    with tab_umaps:
        if 'X_umap' not in adata_display.obsm:
            st.warning("UMAP 2D no calculado o falló.")
        else:
            st.subheader("UMAP 2D coloreado por Clústeres de Leiden")
            if 'leiden_clusters' in adata_display.obs:
                fig_umap_c, ax_c = plt.subplots()
                sc.pl.umap(adata_display, color='leiden_clusters', legend_loc='on data', ax=ax_c, show=False, 
                           title=f"Res: {st.session_state.leiden_res}", size=st.session_state.plot_point_size, palette=st.session_state.plot_palette)
                st.pyplot(fig_umap_c)
                st.download_button("UMAP 2D Clústeres (PNG)", fig_to_bytes(fig_umap_c), "umap_clusters.png", key="dl_umc_v3")
                plt.close(fig_umap_c)

            if 'sample' in adata_display.obs:
                st.subheader("UMAP 2D coloreado por Muestra")
                fig_umap_s, ax_s = plt.subplots()
                sc.pl.umap(adata_display, color='sample', ax=ax_s, show=False, title="Por Muestra", size=st.session_state.plot_point_size, palette=st.session_state.plot_palette)
                st.pyplot(fig_umap_s)
                st.download_button("UMAP 2D Muestra (PNG)", fig_to_bytes(fig_umap_s), "umap_sample.png", key="dl_ums_v3")
                plt.close(fig_umap_s)
            
            if st.session_state.calc_umap_3d and 'X_umap_3d' in adata_display.obsm and 'leiden_clusters' in adata_display.obs:
                st.subheader("UMAP 3D Interactivo por Clústeres")
                try:
                    umap_3d_coords = adata_display.obsm['X_umap_3d']
                    df_umap3d = pd.DataFrame({
                        'UMAP1': umap_3d_coords[:, 0], 'UMAP2': umap_3d_coords[:, 1], 'UMAP3': umap_3d_coords[:, 2],
                        'Cluster': adata_display.obs['leiden_clusters'].astype(str),
                        'Muestra': adata_display.obs['sample'].astype(str)
                    })
                    n_clusters_3d = adata_display.obs['leiden_clusters'].nunique()
                    color_discrete_sequence_3d = px.colors.qualitative.Plotly if n_clusters_3d <= 10 else px.colors.qualitative.Alphabet 
                    
                    fig_3d = px.scatter_3d(
                        df_umap3d, x='UMAP1', y='UMAP2', z='UMAP3', color='Cluster',
                        hover_data=['Muestra'], title="UMAP 3D por Clústeres Leiden",
                        color_discrete_sequence=color_discrete_sequence_3d
                    )
                    fig_3d.update_traces(marker=dict(size=max(1, int(st.session_state.plot_point_size / 10))))
                    st.plotly_chart(fig_3d, use_container_width=True)
                except Exception as e_plot3d_tab:
                    st.error(f"Error generando UMAP 3D: {e_plot3d_tab}")
            elif st.session_state.calc_umap_3d:
                 st.info("UMAP 3D fue seleccionado pero no se pudo calcular o los datos necesarios faltan.")

            st.subheader("UMAPs 2D por Muestra (Facetado, Coloreado por Clúster)")
            if 'sample' in adata_display.obs and 'leiden_clusters' in adata_display.obs:
                try:
                    unique_samples_facet = sorted(adata_display.obs['sample'].astype('category').cat.categories.tolist())
                    n_samples_facet = len(unique_samples_facet)
                    if n_samples_facet > 0:
                        cols_facet = min(n_samples_facet, 3); rows_facet = (n_samples_facet + cols_facet - 1) // cols_facet
                        fig_facet, axes_facet = plt.subplots(rows_facet, cols_facet, figsize=(cols_facet * 5.5, rows_facet * 5), squeeze=False) # Ajustar tamaño
                        axes_flat_facet = axes_facet.flatten()
                        idx_facet = 0 
                        for idx_facet, sample_val_facet in enumerate(unique_samples_facet):
                            if idx_facet < len(axes_flat_facet):
                                ax_curr_facet = axes_flat_facet[idx_facet]
                                adata_subset_facet = adata_display[adata_display.obs['sample'] == sample_val_facet].copy()
                                if not adata_subset_facet.obs.empty:
                                    sc.pl.umap(adata_subset_facet, color='leiden_clusters', ax=ax_curr_facet, show=False, 
                                               title=f"Muestra: {sample_val_facet}", legend_loc='on data' if idx_facet == 0 and n_samples_facet > 1 else None, # Leyenda solo si hay espacio y es el primero
                                               legend_fontsize=6, size=st.session_state.plot_point_size, palette=st.session_state.plot_palette)
                                else:
                                    ax_curr_facet.text(0.5, 0.5, f"M: {sample_val_facet}\n(Sin células)", ha='center', va='center'); ax_curr_facet.set_xticks([]); ax_curr_facet.set_yticks([])
                        for j_ax_empty in range(idx_facet + 1, len(axes_flat_facet)): fig_facet.delaxes(axes_flat_facet[j_ax_empty])
                        plt.tight_layout(); st.pyplot(fig_facet)
                        st.download_button("Descargar UMAPs Facetados (PNG)", fig_to_bytes(fig_facet), "umaps_faceted.png", key="dl_umaps_facet_v3")
                        plt.close(fig_facet)
                except Exception as e_facet: st.error(f"Error UMAPs facetados: {e_facet}")

    with tab_markers:
        if st.session_state.marker_genes_df is not None and not st.session_state.marker_genes_df.empty:
            st.subheader(f"Top {st.session_state.n_top_markers} Genes Marcadores por Clúster")
            st.dataframe(st.session_state.marker_genes_df)
            st.download_button(
                "Descargar Marcadores (CSV)", st.session_state.marker_genes_df.to_csv(index=False).encode('utf-8'),
                "cluster_markers.csv", "text/csv", key="dl_markers_v3"
            )
            
            st.subheader("Dot Plot de Genes Marcadores")
            df_markers_for_plot = st.session_state.marker_genes_df
            # Slider para que el usuario elija cuántos genes por cluster para el dotplot de marcadores
            top_n_per_cluster_dotplot_val = st.slider("Nº top marcadores/clúster para Dot Plot", 1, 
                                                      min(5, st.session_state.n_top_markers), # Max 5 o n_top_markers
                                                      max(1,min(2,st.session_state.n_top_markers)), # Default a 2 o 1
                                                      key="dotplot_n_markers_slider_tab")
            genes_for_marker_dotplot_list = []
            if 'Cluster' in df_markers_for_plot.columns and 'Rank' in df_markers_for_plot.columns and 'Gene' in df_markers_for_plot.columns:
                for cluster_id_dotplot in sorted(df_markers_for_plot['Cluster'].astype(str).unique()):
                    cluster_specific_markers = df_markers_for_plot[df_markers_for_plot['Cluster'] == cluster_id_dotplot]
                    genes_for_marker_dotplot_list.extend(cluster_specific_markers.nsmallest(top_n_per_cluster_dotplot_val, 'Rank')['Gene'].tolist())
                unique_genes_for_dotplot = list(dict.fromkeys(genes_for_marker_dotplot_list))

                if unique_genes_for_dotplot and 'leiden_clusters' in adata_display.obs:
                    try:
                        num_clusters_dotplot = adata_display.obs['leiden_clusters'].nunique()
                        fig_dot_markers, ax_dot_markers = plt.subplots(figsize=(max(8, len(unique_genes_for_dotplot) * 0.6), max(5, num_clusters_dotplot * 0.45)))
                        sc.pl.dotplot(adata_display, unique_genes_for_dotplot, groupby='leiden_clusters', ax=ax_dot_markers, show=False, standard_scale='var', use_raw=False)
                        plt.xticks(rotation=90); plt.tight_layout(); st.pyplot(fig_dot_markers)
                        st.download_button("Descargar Dot Plot Marcadores (PNG)", fig_to_bytes(fig_dot_markers), "dotplot_markers.png", key="dl_dotplot_markers_v3")
                        plt.close(fig_dot_markers)
                    except Exception as e_dot_m: st.error(f"Error dot plot marcadores: {e_dot_m}")
                else: st.info("No hay genes marcadores para el dot plot o faltan clusters.")
            else: st.info("Columnas requeridas no encontradas para dot plot de marcadores.")
        else: st.info("No se han calculado genes marcadores.")


    with tab_heatmap: # Nueva pestaña para Heatmap de Marcadores
            st.subheader(f"Heatmap de Top {st.session_state.heatmap_top_n_genes} Genes Marcadores por Clúster")
            if st.session_state.marker_genes_df is not None and not st.session_state.marker_genes_df.empty:
                heatmap_genes_list = []
                df_markers_for_plot_hm = st.session_state.marker_genes_df
                if 'Cluster' in df_markers_for_plot_hm.columns and 'Rank' in df_markers_for_plot_hm.columns and 'Gene' in df_markers_for_plot_hm.columns:
                    for cluster_id_hm in sorted(df_markers_for_plot_hm['Cluster'].astype(str).unique()):
                        cluster_markers_hm = df_markers_for_plot_hm[df_markers_for_plot_hm['Cluster'] == cluster_id_hm]
                        heatmap_genes_list.extend(cluster_markers_hm.nsmallest(st.session_state.heatmap_top_n_genes, 'Rank')['Gene'].tolist())
                    unique_heatmap_genes = list(dict.fromkeys(heatmap_genes_list))
                else:
                    unique_heatmap_genes = []

                if unique_heatmap_genes:
                    try:
                        st.write(f"Generando heatmap para {len(unique_heatmap_genes)} genes únicos.")
                        adata_for_heatmap = adata_display.copy() 
                        genes_present_in_adata_for_heatmap = [g for g in unique_heatmap_genes if g in adata_for_heatmap.var_names]
                        
                        if not genes_present_in_adata_for_heatmap:
                            st.warning("Ninguno de los genes marcadores seleccionados para el heatmap se encontró en los datos.")
                        else:
                            # --- Calcular dendrograma explícitamente con la key por defecto ---
                            dendrogram_calculated_successfully = False
                            default_dendrogram_key = f"dendrogram_{'leiden_clusters'}"

                            if 'X_pca_hvg' in adata_for_heatmap.obsm:
                                print(f"DEBUG Heatmap: Usando 'X_pca_hvg' para el dendrograma, guardando en '{default_dendrogram_key}'.")
                                sc.tl.dendrogram(
                                    adata_for_heatmap, 
                                    groupby='leiden_clusters', 
                                    use_rep='X_pca_hvg', 
                                    key_added=default_dendrogram_key
                                )
                                dendrogram_calculated_successfully = default_dendrogram_key in adata_for_heatmap.uns
                            elif 'X_pca' in adata_for_heatmap.obsm : 
                                print(f"DEBUG Heatmap: Usando 'X_pca' para el dendrograma, guardando en '{default_dendrogram_key}'.")
                                sc.tl.dendrogram(
                                    adata_for_heatmap, 
                                    groupby='leiden_clusters', 
                                    use_rep='X_pca', 
                                    key_added=default_dendrogram_key
                                )
                                dendrogram_calculated_successfully = default_dendrogram_key in adata_for_heatmap.uns
                            else:
                                print("DEBUG Heatmap: No se encontró X_pca_hvg o X_pca. No se calculará dendrograma explícito para sc.pl.heatmap.")
                            
                            print(f"DEBUG Heatmap: dendrogram_calculated_successfully = {dendrogram_calculated_successfully}")

                            # --- Lógica de ploteo del Heatmap ---
                            if dendrogram_calculated_successfully:
                                print("DEBUG Heatmap: Dendrograma calculado. Intentando plot SIN pasar 'ax' y guardando en archivo.")
                                
                                # Usar un nombre de archivo base para el guardado temporal
                                # Scanpy le añadirá prefijos y directorios
                                base_temp_filename_for_save = "scanpy_heatmap_temp.png"
                                
                                # Llamar a sc.pl.heatmap sin 'ax' y con 'save'
                                # show=False es crucial
                                sc.pl.heatmap(
                                    adata_for_heatmap, 
                                    genes_present_in_adata_for_heatmap, 
                                    groupby='leiden_clusters', 
                                    cmap=st.session_state.plot_palette, 
                                    standard_scale='var', 
                                    dendrogram=True, # Debería encontrar la clave por defecto
                                    show=False,
                                    save=base_temp_filename_for_save # Pasar solo el sufijo/nombre base
                                )
                                
                                # Scanpy guarda en "figures/heatmap_" + save_arg
                                # Necesitamos construir la ruta esperada
                                expected_figures_dir = "figures"
                                expected_filename = f"heatmap{base_temp_filename_for_save}" # Scanpy añade 'heatmap'
                                expected_filepath = os.path.join(expected_figures_dir, expected_filename)
                                
                                print(f"DEBUG Heatmap: Buscando archivo guardado en '{expected_filepath}'")

                                if os.path.exists(expected_filepath):
                                    st.image(expected_filepath) # Mostrar la imagen guardada
                                    with open(expected_filepath, "rb") as f_img_heatmap:
                                        st.download_button(
                                            "Descargar Heatmap Marcadores (PNG)", 
                                            f_img_heatmap.read(), 
                                            "heatmap_markers.png", # Nombre de descarga para el usuario
                                            "image/png", 
                                            key="dl_heatmap_markers_png_via_save"
                                        )
                                    os.remove(expected_filepath) # Limpiar el archivo temporal
                                else:
                                    st.error(f"No se pudo generar o encontrar el archivo temporal del heatmap en '{expected_filepath}'. Revisa si el directorio 'figures' fue creado y el archivo existe.")
                                    # Listar contenido de 'figures' si existe, para depuración
                                    if os.path.exists(expected_figures_dir):
                                        st.text(f"Contenido de '{expected_figures_dir}': {os.listdir(expected_figures_dir)}")


                            else: # Si no hay dendrograma, podemos usar 'ax' de forma segura
                                st.info("No se pudo calcular el dendrograma (probablemente falta X_pca_hvg). Mostrando heatmap sin dendrograma.")
                                fig_heatmap_no_dendro, ax_heatmap_no_dendro = plt.subplots(figsize=(10, max(6, len(genes_present_in_adata_for_heatmap) * 0.3)))
                                sc.pl.heatmap(
                                    adata_for_heatmap, 
                                    genes_present_in_adata_for_heatmap, 
                                    groupby='leiden_clusters', 
                                    cmap=st.session_state.plot_palette, 
                                    standard_scale='var', 
                                    dendrogram=False, # Forzar a False
                                    ax=ax_heatmap_no_dendro, 
                                    show=False
                                )
                                st.pyplot(fig_heatmap_no_dendro)
                                st.download_button("Descargar Heatmap (sin dendro) (PNG)", fig_to_bytes(fig_heatmap_no_dendro), "heatmap_no_dendro.png", key="dl_heatmap_no_dendro_png_v2")
                                plt.close(fig_heatmap_no_dendro)

                    except Exception as e_heatmap_tab:
                        st.error(f"Error generando heatmap: {e_heatmap_tab}")
                        st.error(traceback.format_exc())
                else:
                    st.info("No hay genes marcadores seleccionados para el heatmap.")
            else:
                st.info("No se han calculado genes marcadores.")
    

    with tab_qc:
        # ... (código de QC como lo tenías, usando adata_display) ...
        st.subheader("Plots de Control de Calidad (Post-Filtrado)")
        if 'sample' in adata_display.obs:
            qc_metrics_list = ['n_genes_by_counts', 'total_counts', 'pct_counts_mt']
            qc_titles_display = ["Nº Genes Detectados por Célula", "Nº Total de Cuentas por Célula", "% Cuentas Mitocondriales"]
            for metric_key, metric_title in zip(qc_metrics_list, qc_titles_display):
                if metric_key in adata_display.obs.columns:
                    st.markdown(f"#### {metric_title} (por Muestra)")
                    try:
                        n_samples_qc = adata_display.obs['sample'].nunique()
                        fig_qc_violin, ax_qc_violin = plt.subplots(figsize=(max(6, n_samples_qc * 1.2), 5))
                        sc.pl.violin(adata_display, keys=metric_key, groupby='sample', rotation=45, ax=ax_qc_violin, show=False, cut=0, use_raw=False)
                        ax_qc_violin.set_xlabel("Muestra"); plt.tight_layout(); st.pyplot(fig_qc_violin)
                        st.download_button(f"Descargar {metric_title.replace(' ', '_')} (PNG)", fig_to_bytes(fig_qc_violin), f"qc_violin_{metric_key}.png", "image/png", key=f"dl_qc_vln_{metric_key}_v3")
                        plt.close(fig_qc_violin)
                    except Exception as e_qc_plot: st.error(f"Error violin QC para {metric_title}: {e_qc_plot}")
                else: st.warning(f"Métrica QC '{metric_key}' no encontrada.")
        else: st.warning("Columna 'sample' no encontrada para plots de QC.")


    with tab_dea:
        # ... (código de DEA como lo tenías, usando adata_display si es necesario para alguna configuración, pero el DEA en sí usa adata_for_dea_preview) ...
        st.subheader("Resultados del Análisis de Expresión Diferencial")
        if st.session_state.dea_results_df is not None and not st.session_state.dea_results_df.empty:
            st.markdown(f"**Comparación Actual:** `{st.session_state.dea_comparison_str}`")
            st.dataframe(st.session_state.dea_results_df.head(st.session_state.dea_n_genes_display))
            
            csv_dea_filename_safe = "".join(c if c.isalnum() or c in (' ', '_', '-') else '_' for c in st.session_state.dea_comparison_str).rstrip()
            csv_dea_filename = f"dea_results_{csv_dea_filename_safe.replace(' vs ','_vs_').replace(' ','_')}.csv"
            st.download_button(
                "Descargar Tabla DEA Completa (CSV)", 
                st.session_state.dea_results_df.to_csv(index=False).encode('utf-8'),
                csv_dea_filename, "text/csv", key="dl_dea_table_csv_v3"
            )
            
            st.markdown("#### Volcano Plot Interactivo")
            try:
                df_volcano_plot = st.session_state.dea_results_df.copy()
                df_volcano_plot['Significancia'] = 'No Significativo'
                up_criteria = (df_volcano_plot['Log2FC'] > st.session_state.dea_lfc_cutoff) & (df_volcano_plot['P-Value Adj'] < st.session_state.dea_pval_cutoff)
                down_criteria = (df_volcano_plot['Log2FC'] < -st.session_state.dea_lfc_cutoff) & (df_volcano_plot['P-Value Adj'] < st.session_state.dea_pval_cutoff)
                df_volcano_plot.loc[up_criteria, 'Significancia'] = 'Upregulated'
                df_volcano_plot.loc[down_criteria, 'Significancia'] = 'Downregulated'
                df_volcano_plot['-log10 P-Value Adj'] = -np.log10(df_volcano_plot['P-Value Adj'] + 1e-300) 
                
                fig_px_volcano = px.scatter(
                    df_volcano_plot, x='Log2FC', y='-log10 P-Value Adj', color='Significancia',
                    color_discrete_map={'Upregulated': 'red', 'Downregulated': 'blue', 'No Significativo': 'grey'},
                    hover_data=['Gene', 'Log2FC', 'P-Value Adj', 'Scores'], 
                    title=f"Volcano Plot: {st.session_state.dea_comparison_str}",
                    labels={'Log2FC': 'Log2 Fold Change', '-log10 P-Value Adj': '-log10(P-valor Ajustado)'}
                )
                fig_px_volcano.add_hline(y=-np.log10(st.session_state.dea_pval_cutoff + 1e-300), line_dash="dash", line_color="black", annotation_text=f"P.adj={st.session_state.dea_pval_cutoff}")
                fig_px_volcano.add_vline(x=st.session_state.dea_lfc_cutoff, line_dash="dash", line_color="black", annotation_text=f"LFC={st.session_state.dea_lfc_cutoff}")
                fig_px_volcano.add_vline(x=-st.session_state.dea_lfc_cutoff, line_dash="dash", line_color="black", annotation_text=f"LFC={-st.session_state.dea_lfc_cutoff}")
                st.plotly_chart(fig_px_volcano, use_container_width=True)
                
                html_volcano_buffer = io.StringIO()
                fig_px_volcano.write_html(html_volcano_buffer)
                html_volcano_filename_safe = "".join(c if c.isalnum() or c in (' ', '_', '-') else '_' for c in st.session_state.dea_comparison_str).rstrip()
                html_volcano_filename = f"volcano_plot_{html_volcano_filename_safe.replace(' vs ','_vs_').replace(' ','_')}.html"
                st.download_button("Descargar Volcano Plot (HTML)", html_volcano_buffer.getvalue(), html_volcano_filename, "text/html", key="dl_volcano_html_v3")
            except Exception as e_volc:
                st.error(f"Error generando Volcano Plot: {e_volc}")
                st.error(traceback.format_exc())
        elif st.session_state.analysis_done:
            st.info("No hay resultados de DEA para mostrar. Ejecuta el Análisis Diferencial.")

    with tab_gene_explorer:
        # ... (código del explorador de genes como lo tenías, usando adata_display) ...
        st.subheader("Visualización de Expresión para Genes Específicos")
        if not genes_to_visualize_list:
            st.info("Ingresa nombres de genes válidos para visualizarlos.")
        else:
            st.write(f"Mostrando expresión para: **{', '.join(genes_to_visualize_list)}**")
            if 'X_umap' not in adata_display.obsm:
                st.warning("Plots UMAP no disponibles ya que UMAP no se calculó o falló.")
            else:
                st.markdown("#### UMAPs coloreados por Expresión Génica")
                n_genes_umap_plot = len(genes_to_visualize_list)
                cols_genes_umap = min(n_genes_umap_plot, 3); rows_genes_umap = (n_genes_umap_plot + cols_genes_umap - 1) // cols_genes_umap
                if n_genes_umap_plot > 0:
                    fig_ge_umaps, axes_ge_umaps = plt.subplots(rows_genes_umap, cols_genes_umap, figsize=(cols_genes_umap * 5, rows_genes_umap * 4.5), squeeze=False)
                    axes_flat_ge_umaps = axes_ge_umaps.flatten()
                    idx_ge_plot = 0 
                    for idx_ge_plot, gene_name_plot in enumerate(genes_to_visualize_list):
                        if idx_ge_plot < len(axes_flat_ge_umaps):
                            ax_ge_curr = axes_flat_ge_umaps[idx_ge_plot]
                            try:
                                sc.pl.umap(adata_display, color=gene_name_plot, ax=ax_ge_curr, show=False, title=gene_name_plot, cmap='viridis', use_raw=False, size=st.session_state.plot_point_size)
                            except Exception as e_ge_umap_plot: 
                                ax_ge_curr.text(0.5, 0.5, f"Error plot\n{gene_name_plot}", ha='center', va='center', color='red'); ax_ge_curr.set_xticks([]); ax_ge_curr.set_yticks([])
                                print(f"Error ploteando UMAP para gen {gene_name_plot}: {e_ge_umap_plot}") 
                    for j_ge_empty_ax in range(idx_ge_plot + 1, len(axes_flat_ge_umaps)): fig_ge_umaps.delaxes(axes_flat_ge_umaps[j_ge_empty_ax])
                    plt.tight_layout(); st.pyplot(fig_ge_umaps)
                    st.download_button("Descargar UMAPs de Genes (PNG)", fig_to_bytes(fig_ge_umaps), "gene_explorer_umaps.png", "image/png", key="dl_ge_umaps_png_v3")
                    plt.close(fig_ge_umaps)
            
            if 'leiden_clusters' in adata_display.obs and genes_to_visualize_list:
                st.markdown("#### Diagramas de Violín por Clúster de Leiden")
                try:
                    genes_for_violin_cluster = genes_to_visualize_list[:min(5, len(genes_to_visualize_list))]
                    n_clusters_violin = adata_display.obs['leiden_clusters'].nunique()
                    fig_ge_violins_cl, ax_ge_violins_cl = plt.subplots(figsize=(max(7, n_clusters_violin * 0.8), 5))
                    sc.pl.violin(adata_display, keys=genes_for_violin_cluster, groupby='leiden_clusters', rotation=45, ax=ax_ge_violins_cl, show=False, use_raw=False, cut=0, multi_panel=len(genes_for_violin_cluster)>1)
                    plt.tight_layout(); st.pyplot(fig_ge_violins_cl)
                    st.download_button("Violines por Clúster (PNG)", fig_to_bytes(fig_ge_violins_cl), "ge_violins_cluster.png", key="dl_ge_vln_cl_v3")
                    plt.close(fig_ge_violins_cl)
                except Exception as e_ge_vln_cl: st.error(f"Error violines por clúster: {e_ge_vln_cl}")
            
            if 'condition_temp_dea' in adata_display.obs and adata_display.obs['condition_temp_dea'].nunique() > 1 and genes_to_visualize_list:
                st.markdown("#### Diagramas de Violín por Condición (definida en DEA)")
                try:
                    genes_for_violin_cond = genes_to_visualize_list[:min(5, len(genes_to_visualize_list))]
                    n_conditions_violin = adata_display.obs['condition_temp_dea'].nunique()
                    fig_ge_violins_co, ax_ge_violins_co = plt.subplots(figsize=(max(7, n_conditions_violin * 1.2), 5))
                    sc.pl.violin(adata_display, keys=genes_for_violin_cond, groupby='condition_temp_dea', rotation=45, ax=ax_ge_violins_co, show=False, use_raw=False, cut=0, multi_panel=len(genes_for_violin_cond)>1)
                    plt.tight_layout(); st.pyplot(fig_ge_violins_co)
                    st.download_button("Violines por Condición (PNG)", fig_to_bytes(fig_ge_violins_co), "ge_violins_condition.png", key="dl_ge_vln_co_v2") # Key diferente
                    plt.close(fig_ge_violins_co)
                except Exception as e_ge_vln_co: st.error(f"Error violines por condición: {e_ge_vln_co}")


            if len(genes_to_visualize_list) > 0 and 'leiden_clusters' in adata_display.obs:
                st.markdown("#### Dot Plot de Genes Seleccionados por Clúster")
                try:
                    n_clusters_dot_ge = adata_display.obs['leiden_clusters'].nunique()
                    fig_ge_dotplot, ax_ge_dotplot = plt.subplots(figsize=(max(8, len(genes_to_visualize_list) * 0.7), max(5, n_clusters_dot_ge * 0.5)))
                    sc.pl.dotplot(adata_display, genes_to_visualize_list, groupby='leiden_clusters', ax=ax_ge_dotplot, show=False, standard_scale='var', use_raw=False)
                    plt.xticks(rotation=90); plt.tight_layout(); st.pyplot(fig_ge_dotplot)
                    st.download_button("Descargar Dot Plot de Genes (PNG)", fig_to_bytes(fig_ge_dotplot), "ge_dotplot.png", key="dl_ge_dotplot_v2") # Key diferente
                    plt.close(fig_ge_dotplot)
                except Exception as e_ge_dot: st.error(f"Error dot plot genes seleccionados: {e_ge_dot}")


    with tab_info:
        # ... (código de info como lo tenías, usando adata_display) ...
        st.subheader("Información del Dataset Procesado")
        st.write(f"Total de Células (post-QC): {adata_display.n_obs}")
        st.write(f"Total de Genes (post-QC): {adata_display.n_vars}")
        if st.session_state.adata_hvg_subset is not None: 
            st.write(f"Número de HVGs usados para PCA/Vecinos/UMAP: {st.session_state.adata_hvg_subset.n_vars}")
        if 'sample' in adata_display.obs:
            st.write("Distribución de células por muestra original:")
            st.dataframe(adata_display.obs['sample'].value_counts())
        if 'leiden_clusters' in adata_display.obs:
            st.write("Distribución de células por clúster de Leiden:")
            st.dataframe(adata_display.obs['leiden_clusters'].value_counts().sort_index())
        st.write("Primeras 5 filas de Metadatos de Células (`.obs`):"); st.dataframe(adata_display.obs.head())
        st.write("Primeras 5 filas de Metadatos de Genes (`.var`):"); st.dataframe(adata_display.var.head())
else:
    if st.session_state.adata_raw is None and st.session_state.num_samples > 0 and not all_files_provided:
        st.info("Bienvenido. Por favor, sube todos los archivos para las muestras especificadas y luego haz clic en 'Cargar y Concatenar Datos'.")
    elif st.session_state.num_samples < 1 : 
        st.info("Bienvenido. Ajusta el 'Número de muestras' en la barra lateral para comenzar.")
    elif st.session_state.adata_raw is None and st.session_state.num_samples > 0 and all_files_provided:
        st.info("Archivos listos. Haz clic en 'Cargar y Concatenar Datos' en la barra lateral.")
    elif st.session_state.adata_raw is not None and not st.session_state.analysis_done:
        st.info("Datos cargados. Haz clic en 'Ejecutar Pipeline Principal' en la barra lateral para iniciar el análisis.")
    else: 
        st.info("Bienvenido al Analizador Interactivo de scRNA-seq. Configura tus muestras y parámetros en la barra lateral izquierda para comenzar.")

st.sidebar.markdown("---")
# Sección para generar reporte (muy básica, solo HTML con parámetros)
if st.session_state.analysis_done and st.session_state.adata_processed:
    if st.sidebar.button("Generar Reporte Básico (HTML)", key="generate_report_btn_sidebar"):
        try:
            report_html_parts = ["<html><head><title>Reporte scRNA-seq</title><style>body{font-family: sans-serif;} ul{list-style-type: none; padding-left: 0;} li{margin-bottom: 5px;}</style></head><body>"]
            report_html_parts.append("<h1>Reporte de Análisis Single-Cell RNA-seq</h1>")
            report_html_parts.append(f"<p>Fecha de Generación: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}</p>")
            
            report_html_parts.append("<h2>Parámetros del Pipeline Clave:</h2><ul>")
            params_to_report_list = ["min_genes", "min_cells", "mito_prefix", "max_mito_pct", "n_top_genes_hvg", 
                                "n_pcs", "n_neighbors_val", "leiden_res", "leiden_flavor", 
                                "umap_n_neighbors", "umap_min_dist", "umap_init_pos", "calc_umap_3d"]
            for p_key_report in params_to_report_list:
                report_html_parts.append(f"<li><b>{p_key_report.replace('_', ' ').title()}:</b> {st.session_state.get(p_key_report, 'N/A')}</li>")
            report_html_parts.append("</ul>")

            adata_report_info = st.session_state.adata_processed
            report_html_parts.append("<h2>Resumen del Dataset Procesado:</h2><ul>")
            report_html_parts.append(f"<li>Células: {adata_report_info.n_obs}</li>")
            report_html_parts.append(f"<li>Genes: {adata_report_info.n_vars}</li>")
            if 'sample' in adata_report_info.obs:
                 report_html_parts.append(f"<li>Muestras: {', '.join(sorted(adata_report_info.obs['sample'].unique()))}</li>")
            if 'leiden_clusters' in adata_report_info.obs:
                 report_html_parts.append(f"<li>Clusters Leiden Encontrados: {adata_report_info.obs['leiden_clusters'].nunique()}</li>")
            report_html_parts.append("</ul>")
            
            # (Idea para el futuro: Incrustar plots como imágenes base64)
            
            report_html_parts.append("</body></html>")
            final_report_html_content = "".join(report_html_parts)
            
            st.sidebar.download_button(
                "Descargar Reporte (HTML)",
                data=final_report_html_content,
                file_name=f"scRNAseq_report_{pd.Timestamp.now().strftime('%Y%m%d_%H%M')}.html",
                mime="text/html",
                key="download_report_html_sidebar_btn"
            )
            st.sidebar.success("Reporte HTML listo para descargar.")
        except Exception as e_report_gen:
            st.sidebar.error(f"Error generando reporte: {e_report_gen}")

# --- Sidebar ---
st.sidebar.header("Configuración de Parámetros")
st.sidebar.subheader("Parámetros de Entrada")
st.sidebar.markdown("### Parámetros de Entrada")
st.sidebar.markdown("Sube los archivos de matriz, características y códigos de barras para cada muestra. Los archivos deben ser válidos y compatibles con 10X Genomics.")
st.sidebar.markdown("#### Archivos de Muestra")
st.sidebar.markdown("Selecciona el número de muestras y sube los archivos correspondientes.")
st.sidebar.markdown("#### Parámetros de Filtrado")
st.sidebar.markdown("Ajusta los parámetros de filtrado para el control de calidad (QC) de las células y genes.")
st.sidebar.markdown("#### Parámetros de Análisis")
st.sidebar.markdown("Ajusta los parámetros para el análisis de datos, incluyendo la selección de genes altamente variables (HVGs), PCA, UMAP y clustering.")
st.sidebar.markdown("#### Parámetros de DEA")
st.sidebar.markdown("Ajusta los parámetros para el análisis diferencial de expresión (DEA).")
st.sidebar.markdown("#### Parámetros de Visualización")
st.sidebar.markdown("Ajusta los parámetros para la visualización de los resultados, incluyendo la selección de genes y la configuración de los plots.")

# Nota sobre dependencias y versiones
st.sidebar.markdown("---")
st.sidebar.info("Analizador scRNA-seq v0.9. Basado en Scanpy.")
st.sidebar.markdown("Si experimentas errores con UMAP (ej: `ValueError: high is out of bounds`), considera usar un entorno con `numpy<2.0` o Python 3.10/3.11.")
st.sidebar.markdown("Se recomienda crear un entorno virtual con versiones compatibles de las bibliotecas (ej: `numpy<2.0` si se experimentan errores con UMAP).")
st.sidebar.markdown("Si tienes problemas, consulta la [documentación de Scanpy](https://scanpy.readthedocs.io/en/stable/) o el [repositorio de GitHub]).")
st.sidebar.markdown("Para más información, visita el [repositorio de GitHub]).")
st.sidebar.markdown("**Desarrollado por:** Pedro Botías - pbotias@ucm.es - https://github.com/pbotiast/scRNASeq")
st.sidebar.markdown("**Licencia:** Licencia MIT - https://opensource.org/licenses/MIT")
st.sidebar.markdown("**Fecha:** [05/05/2025] - [12/05/2025]")  
st.sidebar.markdown("**Versión:** 0.9")
st.sidebar.markdown("**Última Actualización:** 2025-05-12")
st.sidebar.markdown("**Notas:** Esta aplicación es un prototipo y puede contener errores. Usa bajo tu propio riesgo.")
st.sidebar.markdown("**Disclaimer:** Esta aplicación es un prototipo y puede contener errores. Usa bajo tu propio riesgo.")




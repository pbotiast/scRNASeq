# -*- coding: utf-8 -*-
"""
Streamlit app para análisis interactivo de Single-Cell RNA-seq con múltiples muestras.
Versión adaptada para análisis de levadura (*Saccharomyces cerevisiae*)
Incluye carga de SGD, Gene Scoring, Varianza PCA, y mejoras de UX.
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
import gzip
import difflib
from scipy import sparse
import json
import umap  # Librería para UMAP
import traceback
import plotly.express as px
# import plotly.graph_objects as go
# from scipy.stats import zscore
import warnings

# Filtrar warnings específicos
warnings.filterwarnings("ignore", category=UserWarning, module="scanpy.preprocessing._highly_variable_genes")
warnings.filterwarnings("ignore", category=FutureWarning, module="scanpy.plotting._anndata")
warnings.filterwarnings("ignore", category=UserWarning, module="umap.umap_")
warnings.filterwarnings("ignore", category=FutureWarning, module="anndata._core.anndata")


# Configuración de la página
st.set_page_config(layout="wide", initial_sidebar_state="expanded")
st.title("Analizador Interactivo de scRNA-seq - Levadura")

# --- Funciones Auxiliares ---
def validate_10x_files(matrix_file, features_file, barcodes_file):
    files_to_check = {
        "Matrix": (matrix_file, "%%MatrixMarket"),
        "Features": (features_file, None),
        "Barcodes": (barcodes_file, None)
    }
    all_valid = True
    error_messages = []
    for name, (file_obj, expected_start) in files_to_check.items():
        try:
            if file_obj is None:
                raise ValueError(f"Archivo {name} no proporcionado.")
            file_obj.seek(0)
            content_buffer = io.BytesIO(file_obj.getbuffer()) # Leer el buffer una vez

            processed_content_buffer = io.BytesIO()
            if file_obj.name.endswith(".gz"):
                try:
                    decompressed_content = gzip.decompress(content_buffer.read())
                    processed_content_buffer.write(decompressed_content)
                except gzip.BadGzipFile:
                    raise ValueError(f"Archivo {name} ({file_obj.name}) parece .gz pero no se pudo descomprimir.")
            else:
                processed_content_buffer.write(content_buffer.read())
            
            processed_content_buffer.seek(0) # Rebobinar para leer

            if expected_start:
                first_line = processed_content_buffer.readline().decode(errors='ignore')
                if not first_line.startswith(expected_start):
                    raise ValueError(f"({file_obj.name}) no es válido. Esperaba '{expected_start}', obtuvo '{first_line[:50]}...'")
            else: # Validación genérica para TSV
                # Leer algunas líneas para ver si es un TSV válido
                pd.read_csv(processed_content_buffer, sep="\t", header=None, nrows=5, comment='%', on_bad_lines='warn')
            file_obj.seek(0) # Reposicionar el buffer original al final
        except ValueError as ve:
            error_messages.append(f"Error de validación en {name}: {ve}")
            all_valid = False
        except Exception as e:
            error_messages.append(f"Error procesando {name} para validación: {e}")
            # st.error(traceback.format_exc()) # Descomentar para depuración más profunda
            all_valid = False
    
    if not all_valid:
        for msg in error_messages:
            st.error(msg) # Mostrar errores en Streamlit
    return all_valid


def load_10x_data(matrix_file, features_file, barcodes_file, sample_name):
    with tempfile.TemporaryDirectory() as temp_dir:
        file_map = {
            "matrix.mtx": matrix_file,
            "features.tsv": features_file,
            "barcodes.tsv": barcodes_file
        }
        for base_name, uploaded_file in file_map.items():
            filename_in_temp = base_name
            uploaded_file.seek(0)
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

def suggest_genes_yeast(input_genes, adata_var_names_lower_map, sgd_aliases_lower_map=None):
    suggestions = {}
    use_sgd = sgd_aliases_lower_map and len(sgd_aliases_lower_map) > 0
    
    valid_adata_genes_lower = list(adata_var_names_lower_map.keys())

    if use_sgd:
        # print("DEBUG suggest_genes: Usando SGD aliases")
        # Para la búsqueda de coincidencias, usamos los alias de SGD
        valid_genes_search_pool_lower = list(sgd_aliases_lower_map.keys())
        # Para mostrar el nombre estándar de SGD en la sugerencia
        def get_standard_name_for_suggestion(alias_lower):
            return sgd_aliases_lower_map.get(alias_lower, alias_lower.upper())
    else:
        # print("DEBUG suggest_genes: Usando var_names de AnnData")
        valid_genes_search_pool_lower = valid_adata_genes_lower
        def get_standard_name_for_suggestion(var_name_lower):
            return adata_var_names_lower_map.get(var_name_lower, var_name_lower.upper())

    for gene_input_raw in input_genes:
        gene_input_stripped = gene_input_raw.strip()
        if not gene_input_stripped: continue
        gene_input_lower = gene_input_stripped.lower()

        # Si el gen ya está en el AnnData con su nombre original (insensible a mayúsculas), no sugerir.
        if gene_input_lower in adata_var_names_lower_map:
            continue 
        
        # Buscar coincidencias cercanas en el pool de búsqueda (SGD o var_names)
        matches_lower = difflib.get_close_matches(gene_input_lower, valid_genes_search_pool_lower, n=3, cutoff=0.6)
        
        if matches_lower:
            # Obtener los nombres estándar/originales de las coincidencias
            suggested_standard_names = [get_standard_name_for_suggestion(m_lower) for m_lower in matches_lower]
            # Filtrar para asegurar que los nombres sugeridos realmente existan en el AnnData actual para ser ploteables
            final_suggestions_for_adata = [s_name for s_name in suggested_standard_names if s_name.lower() in adata_var_names_lower_map]
            
            if final_suggestions_for_adata:
                 suggestions[gene_input_raw] = list(dict.fromkeys(final_suggestions_for_adata)) # Únicos
    return suggestions

# --- Inicialización de st.session_state ---
default_values = {
    "num_samples": 1, "sample_files": {}, "min_genes": 200, "min_cells": 3, "mito_prefix": "Q0",
    "max_mito_pct": 15, "n_top_genes_hvg": 500, "n_pcs": 15, "n_neighbors_val": 15, "leiden_res": 0.5,
    "n_top_markers": 10, "adata_raw": None, "adata_processed": None, "adata_hvg_subset": None,
    "analysis_done": False, "marker_genes_df": None, "condition_assignments": {},
    "dea_group1": None, "dea_group2": None, "dea_cluster_scope": "Todos los Clústeres",
    "dea_n_genes_display": 25, "dea_lfc_cutoff": 0.5, "dea_pval_cutoff": 0.05,
    "dea_results_df": None, "dea_comparison_str": "", 
    "gene_explorer_input": "",                  
    "gene_score_list_input": "Nombre_Firma_Ejemplo: YAL001C, YAL002W, YAL003W\n# Otra_Firma: YAR001W, YAR002C", # <--- CAMBIADA AQUÍ    
    "gene_score_user_lists_input": "Nombre_Firma_Ejemplo: YAL001C, YAL002W, YAL003W\n# Otra_Firma: YAR001W, YAR002C",                
    "gene_score_name_input": "Custom_Gene_Score",  # <--- ASEGÚRATE QUE ESTA CLAVE ESTÉ EXACTAMENTE ASÍ            
    "gene_scores_calculated": {},       
    "leiden_flavor": "leidenalg",                  
    "umap_init_pos": "random",                  
    "umap_n_neighbors": 15,                     
    "umap_min_dist": 0.5,                       
    "calc_umap_3d": False,                      
    "plot_palette": "tab20",                 
    "plot_point_size": 30,                    
    "heatmap_top_n_genes": 3,                   
    "show_pca_variance": True,
    "sgd_genes_loaded": False, 
    "sgd_gene_aliases_lower_map": {} 
}
for key, value in default_values.items():
    if key not in st.session_state:
        st.session_state[key] = value

# --- Sidebar ---
with st.sidebar:
    st.image("https://raw.githubusercontent.com/pbotiast/scRNASeq/refs/heads/main/scanpy_logo.png", width=150, use_container_width=True) 
    st.header("Configuración del Análisis")
    
    with st.expander("Guardar/Cargar Configuración", expanded=False):
        if st.button("Guardar Configuración Actual", key="save_params_btn_sidebar_v10"): 
            params_to_save = {
                k: st.session_state[k] for k in default_values 
                if not k.startswith("adata") and \
                   k not in ["sample_files", "analysis_done", "marker_genes_df", 
                              "dea_results_df", "gene_scores_calculated", 
                              "sgd_genes_loaded", "sgd_gene_aliases_lower_map"]
            }
            params_json = json.dumps(params_to_save, indent=4)
            st.download_button(
                label="Descargar parámetros (.json)", data=params_json,
                file_name=f"scRNAseq_app_params_{pd.Timestamp.now().strftime('%Y%m%d_%H%M')}.json",
                mime="application/json", key="download_params_json_sidebar_btn_v10" 
            )
        uploaded_params_file_sidebar = st.file_uploader("Cargar Configuración (.json)", type="json", key="upload_params_sidebar_v10") 
        if uploaded_params_file_sidebar is not None:
            try:
                loaded_params = json.load(uploaded_params_file_sidebar)
                for key, value_loaded in loaded_params.items():
                    if key in st.session_state and not key.startswith("adata") and key != "sample_files":
                        st.session_state[key] = value_loaded
                st.success("Parámetros cargados.")
                st.info("Revisa los parámetros y vuelve a ejecutar el pipeline si es necesario.")
            except Exception as e_load_params_sidebar: st.error(f"Error al cargar parámetros: {e_load_params_sidebar}")

    with st.expander("1. Carga de Datos", expanded=True):
        st.session_state.num_samples = st.number_input(
            "Número de muestras", min_value=1, max_value=10, value=st.session_state.num_samples, step=1, key="num_samples_main_v10"
        )
        for i in range(st.session_state.num_samples):
            sample_widget_key_prefix = f"sample_{i}_v10" 
            st.subheader(f"Muestra {i+1}")
            s_name_key = f"sample_name_{i}"
            if s_name_key not in st.session_state.sample_files: st.session_state.sample_files[s_name_key] = f"Muestra{i+1}"
            st.session_state.sample_files[s_name_key] = st.text_input(f"Nombre Muestra {i+1}", st.session_state.sample_files[s_name_key], key=f"{sample_widget_key_prefix}_name")
            for file_type_key, label in zip(["matrix_file", "features_file", "barcodes_file"],["Matrix (.mtx/.mtx.gz)", "Features (.tsv/.tsv.gz)", "Barcodes (.tsv/.tsv.gz)"]):
                full_key = f"{file_type_key}_{i}"
                if full_key not in st.session_state.sample_files: st.session_state.sample_files[full_key] = None
                st.session_state.sample_files[full_key] = st.file_uploader(f"{label}", type=["mtx", "tsv", "gz"], key=f"{sample_widget_key_prefix}_{file_type_key}")
            if i < st.session_state.num_samples - 1: st.markdown("---")

    with st.expander("2. Parámetros del Pipeline"):
        param_tabs = st.tabs(["Filtrado QC", "HVGs & PCA", "UMAP", "Clustering"])
        with param_tabs[0]: 
            st.session_state.min_genes = st.slider("Mín genes/célula", 50, 1000, st.session_state.min_genes, key="ming_final_qc_v10", help="Filtra células con menos de N genes.")
            st.session_state.min_cells = st.slider("Mín células/gen", 1, 50, st.session_state.min_cells, key="minc_final_qc_v10", help="Filtra genes en menos de N células.")
            st.session_state.mito_prefix = st.text_input("Prefijo mitocondrial", st.session_state.mito_prefix, key="mitop_final_qc_v10", help="Ej: 'Q0' o 'QCR' para levadura.")
            st.session_state.max_mito_pct = st.slider("Máx % mito", 1, 100, st.session_state.max_mito_pct, key="maxmito_final_qc_v10", help="Filtra células con alto % mitocondrial.")
        with param_tabs[1]: 
            st.session_state.n_top_genes_hvg = st.slider("Nº HVGs", 50, 2000, st.session_state.n_top_genes_hvg, key="nhvg_final_pca_v10", help="Nº de Genes Altamente Variables.")
            st.session_state.n_pcs = st.slider("Nº PCs", 5, 50, st.session_state.n_pcs, key="npcs_final_pca_v10", help="En levadura, 10–20 PCs suelen ser suficientes. Usado para PCA y Vecinos.")
            st.session_state.n_neighbors_val = st.slider("Nº Vecinos (grafo KNN)", 2, 50, st.session_state.n_neighbors_val, key="nneigh_final_pca_v10", help="Para UMAP y Leiden.")
        with param_tabs[2]: 
            st.session_state.calc_umap_3d = st.checkbox("Calcular UMAP 3D", st.session_state.calc_umap_3d, key="calc3d_final_umap_v10")
            uinit_list_sb_umap = ["spectral", "random", "pca"]
            st.session_state.umap_init_pos = st.selectbox("Inicialización UMAP", uinit_list_sb_umap, 
                                                          index=uinit_list_sb_umap.index(st.session_state.umap_init_pos) if st.session_state.umap_init_pos in uinit_list_sb_umap else 1, # Default a 'random'
                                                          key="uinit_final_umap_v10", help="'random' más estable con NumPy 2.x.")
            st.session_state.umap_n_neighbors = st.slider("Nº Vecinos UMAP (embedding)", 2, 200, st.session_state.umap_n_neighbors, key="uneigh_final_umap_v10")
            st.session_state.umap_min_dist = st.slider("Distancia Mínima UMAP", 0.0, 1.0, st.session_state.umap_min_dist, 0.01, key="umind_final_umap_v10")
        with param_tabs[3]: 
            lflav_list_sb_clu = ["igraph", "leidenalg"]
            st.session_state.leiden_flavor = st.selectbox("Backend Leiden", lflav_list_sb_clu, 
                                                          index=lflav_list_sb_clu.index(st.session_state.leiden_flavor) if st.session_state.leiden_flavor in lflav_list_sb_clu else 0, 
                                                          key="lflav_final_clu_v10", help="'igraph' recomendado.")
            st.session_state.leiden_res = st.slider("Resolución Leiden", 0.1, 1.0, st.session_state.leiden_res, 0.1, key="lres_final_clu_v10", help="Para levadura, valores más bajos suelen ser suficientes. Mayor res = más clústeres.")
            st.session_state.n_top_markers = st.slider("Nº marcadores/clúster", 1, 25, st.session_state.n_top_markers, key="nmark_final_clu_v10")

    with st.expander("Carga de Archivo SGD (Opcional para Levadura)", expanded=False):
        sgd_file_upload = st.file_uploader("Sube SGD_features.tab (genes de levadura)", type=["tab", "tsv", "gz"], key="sgd_upload_main_v10")
        if sgd_file_upload:
            try:
                if sgd_file_upload.name.endswith(".gz"):
                    with gzip.open(sgd_file_upload, 'rt', encoding='utf-8') as f_sgd:
                        sgd_df_upload = pd.read_csv(f_sgd, sep='\t', na_filter=False, comment="#") # Añadir comment='#'
                else:
                    with io.TextIOWrapper(sgd_file_upload, encoding='utf-8') as f_sgd_text: # Para archivos no comprimidos
                        sgd_df_upload = pd.read_csv(f_sgd_text, sep='\t', na_filter=False, comment="#")

                primary_name_col_options = ['Standard gene name', 'Gene name', 'Standard_gene_name', 'Standard gene name']
                systematic_name_col_options = ['Systematic name', 'Feature name', 'Systematic_name', 'Feature name']
                alias_col_options = ['Alias', 'Aliases']
                feature_type_col_sgd_options = ['Feature type', 'Feature_type', 'type']

                def find_col(df_cols, options):
                    for opt in options:
                        if opt in df_cols: return opt
                    return None

                primary_name_col = find_col(sgd_df_upload.columns, primary_name_col_options)
                systematic_name_col = find_col(sgd_df_upload.columns, systematic_name_col_options)
                alias_col = find_col(sgd_df_upload.columns, alias_col_options)
                feature_type_col_sgd = find_col(sgd_df_upload.columns, feature_type_col_sgd_options)

                required_cols_sgd = [primary_name_col, systematic_name_col, feature_type_col_sgd] # Alias es opcional
                if not all(required_cols_sgd):
                    st.error(f"Columnas esenciales no encontradas en el archivo SGD. Se necesitan columnas para nombre estándar/sistemático y tipo de característica. Columnas encontradas: {sgd_df_upload.columns.tolist()}")
                else:
                    sgd_df_filtered = sgd_df_upload[sgd_df_upload[feature_type_col_sgd] == 'ORF'].copy()
                    temp_sgd_map = {}
                    for _, row in sgd_df_filtered.iterrows():
                        standard_name = str(row[primary_name_col]).strip().upper() if pd.notna(row[primary_name_col]) and str(row[primary_name_col]).strip() else ""
                        systematic_name = str(row[systematic_name_col]).strip().upper() if pd.notna(row[systematic_name_col]) and str(row[systematic_name_col]).strip() else ""
                        main_sgd_name_for_map = standard_name if standard_name else systematic_name
                        if not main_sgd_name_for_map: continue
                        
                        aliases_to_add = set([main_sgd_name_for_map])
                        if systematic_name and systematic_name != main_sgd_name_for_map: aliases_to_add.add(systematic_name)
                        if alias_col and pd.notna(row[alias_col]) and str(row[alias_col]).strip():
                            for alias_item in str(row[alias_col]).split('|'):
                                cleaned_alias = alias_item.strip().upper()
                                if cleaned_alias: aliases_to_add.add(cleaned_alias)
                        for alias_final in aliases_to_add:
                            if alias_final: temp_sgd_map[alias_final.lower()] = main_sgd_name_for_map
                    
                    st.session_state.sgd_gene_aliases_lower_map = temp_sgd_map
                    st.session_state.sgd_genes_loaded = True
                    st.success(f"Archivo SGD procesado: {len(st.session_state.sgd_gene_aliases_lower_map)} alias mapeados.")
                    print(f"DEBUG: SGD Map size: {len(st.session_state.sgd_gene_aliases_lower_map)}")
            except Exception as e_sgd_upload:
                st.error(f"Error al procesar archivo SGD: {e_sgd_upload}"); st.error(traceback.format_exc())
                st.session_state.sgd_genes_loaded = False; st.session_state.sgd_gene_aliases_lower_map = {}

    with st.expander("Personalización de Plots", expanded=False):
        pal_list_sb_custom = ["tab10", "tab20", "Set3", "Paired", "viridis", "plasma", "magma", "cividis", "default"]
        st.session_state.plot_palette = st.selectbox("Paleta Colores UMAP/Heatmap", pal_list_sb_custom, 
                                                     index=pal_list_sb_custom.index(st.session_state.plot_palette) if st.session_state.plot_palette in pal_list_sb_custom else 0,
                                                     key="pal_final_custom_v10")
        st.session_state.plot_point_size = st.slider("Tamaño Puntos UMAP 2D", 10, 150, st.session_state.plot_point_size, 5, key="psize_final_custom_v10")
        st.session_state.heatmap_top_n_genes = st.slider("Nº genes/clúster para Heatmap", 1, 10, st.session_state.heatmap_top_n_genes, key="hmn_final_custom_v10")
        st.session_state.show_pca_variance = st.checkbox("Mostrar Varianza PCA en Info", st.session_state.show_pca_variance, key="showpcavar_final_custom_v10")

    if st.session_state.analysis_done and st.session_state.adata_processed is not None:
        with st.expander("3. Análisis Diferencial (DEA)"):
            adata_for_dea_config_sb = st.session_state.adata_processed
            samples_in_adata_sb = sorted(adata_for_dea_config_sb.obs['sample'].unique().tolist())
            st.subheader("Asignar Condiciones")
            
            current_assignments_dea_sb = {s: st.session_state.condition_assignments.get(s, f"Cond_{s.replace(' ','_')}") for s in samples_in_adata_sb}
            for sample_name_dea_config_sb in samples_in_adata_sb:
                current_assignments_dea_sb[sample_name_dea_config_sb] = st.text_input(
                    f"Condición para {sample_name_dea_config_sb}", 
                    value=current_assignments_dea_sb[sample_name_dea_config_sb], 
                    key=f"cond_assign_{sample_name_dea_config_sb}_final_v10"
                )
            st.session_state.condition_assignments = current_assignments_dea_sb

            unique_defined_conditions_sb = sorted(list(set(c for c in st.session_state.condition_assignments.values() if c and c.strip())))

            if len(unique_defined_conditions_sb) >= 2:
                st.subheader("Seleccionar Grupos para Comparación")
                col1_dea_sb, col2_dea_sb = st.columns(2)
                with col1_dea_sb:
                    st.session_state.dea_group1 = st.selectbox("Grupo 1 (Referencia)", unique_defined_conditions_sb, 
                                                               index=unique_defined_conditions_sb.index(st.session_state.dea_group1) if st.session_state.dea_group1 in unique_defined_conditions_sb else 0,
                                                               key="dea_g1_select_final_v10")
                with col2_dea_sb:
                    options_g2_sb = [c for c in unique_defined_conditions_sb if c != st.session_state.dea_group1]
                    if not options_g2_sb and unique_defined_conditions_sb:
                        options_g2_sb = [c for c in unique_defined_conditions_sb if c != st.session_state.dea_group1] or (unique_defined_conditions_sb[1:] if len(unique_defined_conditions_sb)>1 else unique_defined_conditions_sb)
                    st.session_state.dea_group2 = st.selectbox("Grupo 2 (Comparación)", options_g2_sb, 
                                                               index=options_g2_sb.index(st.session_state.dea_group2) if st.session_state.dea_group2 in options_g2_sb and options_g2_sb else 0,
                                                               key="dea_g2_select_final_v10")

                clusters_for_dea_sb = ["Todos los Clústeres"] + sorted(adata_for_dea_config_sb.obs['leiden_clusters'].astype(str).unique().tolist())
                st.session_state.dea_cluster_scope = st.selectbox("Ámbito DEA", clusters_for_dea_sb, 
                                                                  index=clusters_for_dea_sb.index(st.session_state.dea_cluster_scope) if st.session_state.dea_cluster_scope in clusters_for_dea_sb else 0,
                                                                  key="dea_scope_select_final_v10")
                st.session_state.dea_n_genes_display = st.slider("Nº genes DEA a mostrar", 10, 200, st.session_state.dea_n_genes_display, key="dea_ngenes_slider_final_v10")
                st.session_state.dea_lfc_cutoff = st.number_input("Log2FC cutoff (Volcano)", 0.0, value=st.session_state.dea_lfc_cutoff, step=0.1, key="dea_lfc_input_final_v10")
                st.session_state.dea_pval_cutoff = st.number_input("P-adj cutoff (Volcano)", 0.0, 1.0, value=st.session_state.dea_pval_cutoff, step=0.01, format="%.3f", key="dea_pval_input_final_v10")

                adata_for_dea_preview_sb = adata_for_dea_config_sb.copy()
                adata_for_dea_preview_sb.obs['condition_temp_dea'] = adata_for_dea_preview_sb.obs['sample'].map(st.session_state.condition_assignments)
                valid_cells_for_preview_sb = adata_for_dea_preview_sb.obs['condition_temp_dea'].notna()
                if valid_cells_for_preview_sb.any():
                    counts_df_sb = adata_for_dea_preview_sb[valid_cells_for_preview_sb].obs.groupby(['condition_temp_dea', 'leiden_clusters'], observed=True).size().unstack(fill_value=0)
                    st.write("Conteos Células por Condición/Clúster:"); st.dataframe(counts_df_sb)
                else: st.warning("No hay células con condiciones asignadas.")

                if st.button("Ejecutar DEA", key="run_dea_btn_final_v10"):
                    if not st.session_state.dea_group1 or not st.session_state.dea_group2:
                        st.error("Por favor, selecciona ambos grupos para la comparación.")
                    elif st.session_state.dea_group1 == st.session_state.dea_group2:
                        st.error("Los grupos de comparación deben ser diferentes.")
                    else:
                        with st.spinner("Ejecutando DEA..."):
                            try:
                                adata_filtered_for_dea = adata_for_dea_preview_sb[
                                    adata_for_dea_preview_sb.obs['condition_temp_dea'].isin([st.session_state.dea_group1, st.session_state.dea_group2]) &
                                    adata_for_dea_preview_sb.obs['condition_temp_dea'].notna() 
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
    
    all_files_provided_sb_load_btn = True 
    if st.session_state.num_samples > 0:
        for i in range(st.session_state.num_samples):
            if not (st.session_state.sample_files.get(f"matrix_file_{i}") and \
                    st.session_state.sample_files.get(f"features_file_{i}") and \
                    st.session_state.sample_files.get(f"barcodes_file_{i}")):
                all_files_provided_sb_load_btn = False; break
    else: all_files_provided_sb_load_btn = False

    if all_files_provided_sb_load_btn:
        if st.button("Cargar y Concatenar Datos", key="load_concat_btn_main_v10", type="primary"):
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
                        for key_to_reset in ["marker_genes_df", "dea_results_df", "dea_comparison_str", "dea_group1", "dea_group2", "dea_cluster_scope", "gene_scores_calculated"]:
                            st.session_state[key_to_reset] = default_values.get(key_to_reset)

                        st.success(f"Carga completada: {st.session_state.adata_raw.n_obs} células, {st.session_state.adata_raw.n_vars} genes.")
                        print("DEBUG: Muestras en adata_raw.obs['sample']:", st.session_state.adata_raw.obs['sample'].unique())
                    except Exception as e_load:
                        st.error(f"Error durante la carga: {e_load}")
                        st.error(traceback.format_exc())
                        st.session_state.adata_raw = None
            else: # all_valid_globally es False
                st.session_state.adata_raw = None # Asegurar que no se proceda con datos inválidos
                
    elif st.session_state.num_samples > 0 :
        st.warning("Sube todos los archivos para habilitar la carga.")

    # DEBUG: Mostrar estado de adata_raw antes del botón del pipeline
    if "adata_raw" in st.session_state:
        st.sidebar.caption(f"Estado interno adata_raw: {'Existe' if st.session_state.adata_raw is not None else 'No existe / None'}")
    else:
        st.sidebar.caption("Estado interno adata_raw: Clave no inicializada")


    if st.session_state.adata_raw is not None:
        if st.button("Ejecutar Pipeline Principal", key="run_pipeline_btn_main_v10", type="primary"):
            adata = st.session_state.adata_raw.copy() 
            st.session_state.adata_hvg_subset = None 
            st.session_state.adata_processed = None # Resetear antes de empezar
            st.session_state.analysis_done = False
            
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
                    if 'sample' not in adata.obs.columns: raise KeyError("Columna 'sample' no encontrada en adata.obs para HVG batch_key.")
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
                    
                    try: # UMAP 2D
                        reducer_2d = umap.UMAP(
                            n_neighbors=int(st.session_state.umap_n_neighbors), n_components=2, metric='euclidean',
                            min_dist=float(st.session_state.umap_min_dist), init=st.session_state.umap_init_pos,
                            random_state=42, verbose=False 
                        )
                        embedding_2d = reducer_2d.fit_transform(X_pca_for_umap)
                        st.session_state.adata_hvg_subset.obsm['X_umap'] = embedding_2d
                        print("DEBUG UMAP 2D: API directa completada.")
                    except Exception as e_umap_direct_api_2d: # Renombrar para evitar conflicto
                        st.warning(f"UMAP 2D con API directa falló ({e_umap_direct_api_2d}). Intentando con sc.tl.umap...")
                        try:
                            umap_init_fallback_2d = st.session_state.umap_init_pos if st.session_state.umap_init_pos != 'pca' else 'spectral'
                            sc.tl.umap(st.session_state.adata_hvg_subset, 
                                       n_neighbors=int(st.session_state.umap_n_neighbors), 
                                       min_dist=float(st.session_state.umap_min_dist),
                                       init_pos=umap_init_fallback_2d, 
                                       random_state=0) 
                            print("DEBUG UMAP 2D: sc.tl.umap (fallback) completado.")
                        except Exception as e_sc_umap_fallback_2d:
                             st.error(f"Fallback sc.tl.umap para UMAP 2D también falló: {e_sc_umap_fallback_2d}.")
                             st.session_state.adata_hvg_subset.obsm.pop('X_umap', None) 

                    if st.session_state.calc_umap_3d:
                        status_text.text("Calculando UMAP 3D...")
                        try:
                            reducer_3d = umap.UMAP(
                                n_neighbors=int(st.session_state.umap_n_neighbors), n_components=3, metric='euclidean',
                                min_dist=float(st.session_state.umap_min_dist), init=st.session_state.umap_init_pos,
                                random_state=42, verbose=False
                            )
                            embedding_3d = reducer_3d.fit_transform(X_pca_for_umap)
                            st.session_state.adata_hvg_subset.obsm['X_umap_3d'] = embedding_3d
                            print("DEBUG UMAP 3D: API directa completada.")
                        except Exception as e_umap3d_direct_api: # Renombrar
                            st.warning(f"UMAP 3D con API directa falló ({e_umap3d_direct_api}). No se generará UMAP 3D.")
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
                    if 'X_umap_3d' in st.session_state.adata_hvg_subset.obsm: 
                        adata.obsm['X_umap_3d'] = st.session_state.adata_hvg_subset.obsm['X_umap_3d']
                    if 'X_pca' in st.session_state.adata_hvg_subset.obsm: 
                         adata.obsm['X_pca_hvg'] = st.session_state.adata_hvg_subset.obsm['X_pca'] 

                    # Marcadores
                    adata.obs['leiden_clusters'] = adata.obs['leiden_clusters'].astype('category')
                    sc.tl.rank_genes_groups(adata, 'leiden_clusters', method='wilcoxon', key_added='rank_genes_groups', use_raw=False) # Usar key por defecto 'rank_genes_groups'
                    
                    res_markers = adata.uns.get('rank_genes_groups', {}) # Usar key por defecto
                    marker_data_list = []
                    if res_markers and 'names' in res_markers and hasattr(res_markers['names'], 'dtype') and res_markers['names'].dtype.names is not None:
                        cluster_ids_with_names_markers = res_markers['names'].dtype.names
                        for grp_marker in cluster_ids_with_names_markers:
                            if all(isinstance(res_markers.get(field), np.ndarray) and hasattr(res_markers[field], 'dtype') and grp_marker in res_markers[field].dtype.names for field in ['names', 'scores', 'logfoldchanges', 'pvals_adj']):
                                num_avail_markers = len(res_markers['names'][grp_marker])
                                markers_to_fetch = min(st.session_state.n_top_markers, num_avail_markers)
                                for i_marker in range(markers_to_fetch):
                                    if all(i_marker < len(res_markers[field][grp_marker]) for field in ['names', 'scores', 'logfoldchanges', 'pvals_adj']):
                                        marker_data_list.append({
                                            'Cluster': grp_marker, 'Rank': i_marker + 1, 
                                            'Gene': res_markers['names'][grp_marker][i_marker], 
                                            'Score': res_markers['scores'][grp_marker][i_marker],
                                            'Log2FC': res_markers['logfoldchanges'][grp_marker][i_marker], 
                                            'P-Value Adj': res_markers['pvals_adj'][grp_marker][i_marker]
                                        })
                                    else: print(f"Warn: Discrepancia longitud arrays marcadores para clúster '{grp_marker}', índice {i_marker}")
                            else: print(f"Warn: Faltan campos o el grupo '{grp_marker}' no es un campo válido en todos los arrays de marcadores.")
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
                        # Resetear más estados de session_state a sus valores por defecto
                        for key_to_reset, default_val_to_reset in default_values.items():
                            if not key_to_reset.startswith("adata") and \
                               key_to_reset not in ["sample_files", "num_samples", "analysis_done", 
                                                    "gene_explorer_input", "gene_score_user_lists_input", # No resetear inputs de usuario
                                                    "sgd_genes_loaded", "sgd_gene_aliases_lower_map"]: 
                                st.session_state[key_to_reset] = default_val_to_reset
            # --- FIN DEL PIPELINE PRINCIPAL ---
    # --- FIN SIDEBAR ---


# --- Sección de Resultados ---
st.markdown("---") 
st.header("Resultados del Análisis")

st.subheader("🔬 Explorador de Expresión Génica")
st.session_state.gene_explorer_input = st.text_area(
    "Ingresa nombres de genes (separados por coma, espacio o nueva línea):", 
    value=st.session_state.gene_explorer_input, 
    key="gene_explorer_main_input_final_v2", # Nueva key
    height=100,
    help="Escribe los nombres de los genes que deseas visualizar en UMAPs, violines y dot plots."
)

if st.session_state.analysis_done and st.session_state.adata_processed is not None:
    adata_display = st.session_state.adata_processed 
    # Crear mapa de genes en minúscula a su versión original DEL DATASET ACTUAL para búsqueda
    valid_genes_lower_map_display = {g.lower(): g for g in adata_display.var_names}
    adata_var_names_lower_map = valid_genes_lower_map_display  # Añadido para compatibilidad con el resto del código

    # Botón de descarga de AnnData en la Sidebar
    try:
        with tempfile.NamedTemporaryFile(delete=False, suffix=".h5ad") as tmp_h5ad_file_dl: # Nombre de variable diferente
            filepath_h5ad_dl = tmp_h5ad_file_dl.name
        adata_display.write_h5ad(filepath_h5ad_dl) # Escribir en el path obtenido

        with open(filepath_h5ad_dl, "rb") as f_h5ad_read_dl:
            st.sidebar.download_button( # Moverlo a la sidebar si no estaba ya allí
                "Descargar AnnData Procesado (.h5ad)", 
                f_h5ad_read_dl.read(), 
                f"processed_adata_scRNAseq_{pd.Timestamp.now().strftime('%Y%m%d_%H%M')}.h5ad",
                "application/octet-stream",
                key="download_adata_sidebar_final_v2" # Nueva key
            )
        os.remove(filepath_h5ad_dl) 
    except Exception as e_dl_sidebar_res: 
        st.sidebar.error(f"Error descarga AnnData: {e_dl_sidebar_res}")
        if 'filepath_h5ad_dl' in locals() and os.path.exists(filepath_h5ad_dl):
            try: os.remove(filepath_h5ad_dl)
            except: pass


    # Procesar genes del explorador
    genes_input_list_raw_main = [g.strip() for g in st.session_state.gene_explorer_input.replace(',', ' ').replace('\n', ' ').split() if g.strip()]
    genes_to_visualize_list = [] # Usar este nombre consistentemente
    genes_not_found_in_adata_list_main = []

    for gene_raw_main in genes_input_list_raw_main:
        gene_lower_main = gene_raw_main.lower()
        if gene_lower_main in adata_var_names_lower_map: # Comprobar contra el AnnData actual
            genes_to_visualize_list.append(adata_var_names_lower_map[gene_lower_main])
        else:
            genes_not_found_in_adata_list_main.append(gene_raw_main)
    
    genes_to_visualize_list = list(dict.fromkeys(genes_to_visualize_list)) # Únicos y mantener orden

    if genes_not_found_in_adata_list_main:
        gene_suggestions_main = suggest_genes_yeast(
            genes_not_found_in_adata_list_main, 
            adata_var_names_lower_map, # Mapa de genes del AnnData actual
            st.session_state.sgd_gene_aliases_lower_map if st.session_state.sgd_genes_loaded else None
        )
        warning_msg_main = f"⚠️ Genes no encontrados directamente en los datos: {', '.join(genes_not_found_in_adata_list_main)}."
        if gene_suggestions_main:
            suggestions_str_list_main = [
                f"{gene_orig}: (quizás: {', '.join(sugg_list)}?)" 
                for gene_orig, sugg_list in gene_suggestions_main.items()
            ]
            if suggestions_str_list_main: # Solo añadir si hay sugerencias válidas
                 warning_msg_main += f"\n💡 Sugerencias: {'; '.join(suggestions_str_list_main)}"
        st.warning(warning_msg_main)

    tab_titles_results_main = ["📊 UMAPs", "🔬 Marcadores", "🔥 Heatmap Marcadores", "🎯 Gene Scoring", "🧬 QC", "📈 DEA", "🧬 Explorador Genes", "ℹ️ Info"]
    tabs_results_main_list = st.tabs(tab_titles_results_main) 
    
    tab_umaps_display, tab_markers_display, tab_heatmap_display, tab_gene_scoring_display, \
    tab_qc_display, tab_dea_display, tab_gene_explorer_display, tab_info_display = tabs_results_main_list
    
    # Aplicar paleta de colores seleccionada
    # Es mejor pasar 'palette' directamente a cada función sc.pl para evitar efectos secundarios globales
    # sc.set_figure_params(color_map=st.session_state.plot_palette)

    with tab_umaps_display:
        if 'X_umap' not in adata_display.obsm:
            st.warning("UMAP 2D no calculado o falló. No se pueden mostrar plots UMAP 2D.")
        else:
            st.subheader("UMAP 2D coloreado por Clústeres de Leiden")
            if 'leiden_clusters' in adata_display.obs:
                fig_umap_c_disp, ax_c_disp = plt.subplots(figsize=(7,6)) # figsize
                sc.pl.umap(adata_display, color='leiden_clusters', legend_loc='on data', ax=ax_c_disp, show=False, 
                           title=f"Res: {st.session_state.leiden_res}", size=st.session_state.plot_point_size, 
                           palette=st.session_state.plot_palette if st.session_state.plot_palette != 'default' else None)
                st.pyplot(fig_umap_c_disp)
                st.download_button("UMAP 2D Clústeres (PNG)", fig_to_bytes(fig_umap_c_disp), "umap_clusters.png", key="dl_umc_final")
                plt.close(fig_umap_c_disp)
            else: st.warning("Clusters Leiden no encontrados para plot UMAP.")

            if 'sample' in adata_display.obs:
                st.subheader("UMAP 2D coloreado por Muestra")
                fig_umap_s_disp, ax_s_disp = plt.subplots(figsize=(7,6))
                sc.pl.umap(adata_display, color='sample', ax=ax_s_disp, show=False, title="Por Muestra", 
                           size=st.session_state.plot_point_size, 
                           palette=st.session_state.plot_palette if st.session_state.plot_palette != 'default' else None)
                st.pyplot(fig_umap_s_disp)
                st.download_button("UMAP 2D Muestra (PNG)", fig_to_bytes(fig_umap_s_disp), "umap_sample.png", key="dl_ums_final")
                plt.close(fig_umap_s_disp)
            
            # UMAP 3D
            if st.session_state.calc_umap_3d: # Solo intentar si el usuario lo pidió
                if 'X_umap_3d' in adata_display.obsm and 'leiden_clusters' in adata_display.obs:
                    st.subheader("UMAP 3D Interactivo por Clústeres")
                    try:
                        umap_3d_coords_disp = adata_display.obsm['X_umap_3d']
                        df_umap3d_disp = pd.DataFrame({
                            'UMAP1': umap_3d_coords_disp[:, 0], 'UMAP2': umap_3d_coords_disp[:, 1], 'UMAP3': umap_3d_coords_disp[:, 2],
                            'Cluster': adata_display.obs['leiden_clusters'].astype(str),
                            'Muestra': adata_display.obs['sample'].astype(str) # Añadir muestra a hover
                        })
                        n_clusters_3d_disp = adata_display.obs['leiden_clusters'].nunique()
                        # Usar una paleta de Plotly, o si es 'default', dejar que Plotly elija
                        color_seq_3d = px.colors.qualitative.Plotly if st.session_state.plot_palette == "default" or n_clusters_3d_disp > len(px.colors.qualitative.Plotly) else getattr(px.colors.qualitative, st.session_state.plot_palette, px.colors.qualitative.Plotly)
                        
                        fig_3d_disp = px.scatter_3d(
                            df_umap3d_disp, x='UMAP1', y='UMAP2', z='UMAP3', color='Cluster',
                            hover_data=['Muestra'], title="UMAP 3D por Clústeres Leiden",
                            color_discrete_sequence=color_seq_3d if isinstance(color_seq_3d, list) else None, # Para paletas Plotly
                            color_discrete_map=None # No usar color_map y color_discrete_sequence juntos
                        )
                        fig_3d_disp.update_traces(marker=dict(size=max(1, int(st.session_state.plot_point_size / 15)))) # Tamaño de punto más pequeño para 3D
                        st.plotly_chart(fig_3d_disp, use_container_width=True)
                    except Exception as e_plot3d_tab_disp:
                        st.error(f"Error generando UMAP 3D: {e_plot3d_tab_disp}")
                elif st.session_state.calc_umap_3d: # Si se pidió pero no están los datos
                     st.info("UMAP 3D fue seleccionado en parámetros pero no se pudo calcular o los datos necesarios ('X_umap_3d', 'leiden_clusters') faltan.")

            # UMAPs Facetados
            st.subheader("UMAPs 2D por Muestra (Coloreado por Clúster)")
            if 'sample' in adata_display.obs and 'leiden_clusters' in adata_display.obs:
                try:
                    unique_samples_facet_disp = sorted(adata_display.obs['sample'].astype('category').cat.categories.tolist())
                    n_samples_facet_disp = len(unique_samples_facet_disp)
                    if n_samples_facet_disp > 0:
                        cols_facet_disp = min(n_samples_facet_disp, 3)
                        rows_facet_disp = (n_samples_facet_disp + cols_facet_disp - 1) // cols_facet_disp
                        fig_facet_disp, axes_facet_disp = plt.subplots(rows_facet_disp, cols_facet_disp, 
                                                                    figsize=(cols_facet_disp * 5.5, rows_facet_disp * 5), squeeze=False)
                        axes_flat_facet_disp = axes_facet_disp.flatten()
                        idx_facet_disp = 0 
                        for idx_facet_disp, sample_val_facet_disp in enumerate(unique_samples_facet_disp):
                            if idx_facet_disp < len(axes_flat_facet_disp):
                                ax_curr_facet_disp = axes_flat_facet_disp[idx_facet_disp]
                                adata_subset_facet_disp = adata_display[adata_display.obs['sample'] == sample_val_facet_disp].copy()
                                if not adata_subset_facet_disp.obs.empty and 'X_umap' in adata_subset_facet_disp.obsm: # Verificar X_umap
                                    sc.pl.umap(adata_subset_facet_disp, color='leiden_clusters', ax=ax_curr_facet_disp, show=False, 
                                               title=f"Muestra: {sample_val_facet_disp}", 
                                               legend_loc='on data' if idx_facet_disp == 0 and n_samples_facet_disp > 1 else None, 
                                               legend_fontsize=6, size=st.session_state.plot_point_size, 
                                               palette=st.session_state.plot_palette if st.session_state.plot_palette != 'default' else None)
                                elif not adata_subset_facet_disp.obs.empty:
                                     ax_curr_facet_disp.text(0.5,0.5, f"M: {sample_val_facet_disp}\n(X_umap no disp.)", ha='center',va='center')
                                else:
                                    ax_curr_facet_disp.text(0.5, 0.5, f"M: {sample_val_facet_disp}\n(Sin células)", ha='center', va='center')
                                ax_curr_facet_disp.set_xticks([]); ax_curr_facet_disp.set_yticks([])
                        for j_ax_empty_disp in range(idx_facet_disp + 1, len(axes_flat_facet_disp)): fig_facet_disp.delaxes(axes_flat_facet_disp[j_ax_empty_disp])
                        plt.tight_layout(); st.pyplot(fig_facet_disp)
                        st.download_button("Descargar UMAPs Facetados (PNG)", fig_to_bytes(fig_facet_disp), "umaps_faceted.png", key="dl_umaps_facet_final")
                        plt.close(fig_facet_disp)
                except Exception as e_facet_disp: st.error(f"Error UMAPs facetados: {e_facet_disp}")


    with tab_markers_display:
        st.subheader(f"Top {st.session_state.n_top_markers} Genes Marcadores por Clúster")
        if st.session_state.marker_genes_df is not None and not st.session_state.marker_genes_df.empty:
            
            df_markers_to_display = st.session_state.marker_genes_df.copy() 

            # --- Configuración de Columnas para st.data_editor ---
            column_config_markers = {
                "P-Value": st.column_config.NumberColumn(
                    "P-Valor", # Label opcional para la columna en la UI
                    format="%.2e", # Notación científica con 2 decimales
                    help="P-valor crudo del test de Wilcoxon."
                ),
                "P-Value Adj": st.column_config.NumberColumn(
                    "P-Valor Ajustado",
                    format="%.2e", # Notación científica con 2 decimales
                    help="P-valor ajustado por múltiples comparaciones (Benjamini-Hochberg)."
                ),
                "Score": st.column_config.NumberColumn(
                    "Score",
                    format="%.3f", # 3 decimales flotantes
                ),
                "Log2FC": st.column_config.NumberColumn(
                    "Log2 Fold Change",
                    format="%.3f",
                ),
                "Gene": st.column_config.TextColumn("Gen"), # Asegurar que se trate como texto
                "Cluster": st.column_config.TextColumn("Clúster"),
                "Rank": st.column_config.NumberColumn("Rank", format="%d") # Entero
            }
            
            # Añadir columnas de enlace si el checkbox está activo
            show_links_markers = st.checkbox(
                "Mostrar enlaces de búsqueda para marcadores (SGD, NCBI)", 
                value=True, # Mostrar por defecto
                key="show_marker_links_check_final_v3" 
            )

            if show_links_markers:
                sgd_base_url_marker_tab = "https://www.yeastgenome.org/locus/" # Renombrar para evitar conflicto de scope
                ncbi_base_url_marker_tab = "https://www.ncbi.nlm.nih.gov/gene/?term="
                
                # Asegurarse de que la columna 'Gene' existe antes de intentar aplicar funciones
                if 'Gene' in df_markers_to_display.columns:
                    df_markers_to_display['SGD_Link_Marker'] = df_markers_to_display['Gene'].apply(
                        lambda x: f"{sgd_base_url_marker_tab}{x.upper()}" if pd.notna(x) else ""
                    )
                    df_markers_to_display['NCBI_Link_Marker'] = df_markers_to_display['Gene'].apply(
                        lambda x: f"{ncbi_base_url_marker_tab}{x}%20AND%20saccharomyces%20cerevisiae[orgn]" if pd.notna(x) else ""
                    )
                    
                    column_config_markers["SGD_Link_Marker"] = st.column_config.LinkColumn("SGD", display_text="🔗SGD", width="small")
                    column_config_markers["NCBI_Link_Marker"] = st.column_config.LinkColumn("NCBI", display_text="🔗NCBI", width="small")
                else:
                    st.warning("Columna 'Gene' no encontrada para generar enlaces de búsqueda.")

            st.data_editor(
                df_markers_to_display, 
                height=400, 
                use_container_width=True, 
                num_rows="dynamic",
                column_config=column_config_markers,
                key="marker_genes_data_editor_final_v3" 
            )
            # Descargar el DataFrame original sin las columnas de enlace Markdown si se añadieron solo para display
            st.download_button(
                "Descargar Marcadores (CSV)", 
                st.session_state.marker_genes_df.to_csv(index=False).encode('utf-8'), 
                "cluster_markers.csv", "text/csv", key="dl_markers_final_v2" 
            )
        else:
            st.info("No se han calculado genes marcadores o la tabla está vacía.")

    with tab_heatmap_display:
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





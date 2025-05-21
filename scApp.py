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
    "n_pcs_actually_used_in_pipeline": 30, # O un valor inicial razonable                      
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
                    print(f"DEBUG PCA: Usando n_comps={current_n_pcs_val}")
                    st.session_state.n_pcs_actually_used_in_pipeline = current_n_pcs_val # <--- GUARDAR AQUÍ
                    sc.tl.pca(st.session_state.adata_hvg_subset, svd_solver='arpack', n_comps=current_n_pcs_val, random_state=0)

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
            # --- UMAP 2D POR CLÚSTERES ---
            st.subheader("UMAP 2D por Clústeres de Leiden")
            if 'leiden_clusters' in adata_display.obs:
                
                # --- VERSIÓN INTERACTIVA (PLOTLY) ---
                with st.expander("Ver UMAP 2D Interactivo por Clústeres (Plotly)", expanded=True): # Expandido por defecto
                    try:
                        umap_coords_c_plotly = adata_display.obsm['X_umap']
                        df_umap_c_plotly = pd.DataFrame({
                            'UMAP1': umap_coords_c_plotly[:, 0],
                            'UMAP2': umap_coords_c_plotly[:, 1],
                            'Cluster': adata_display.obs['leiden_clusters'].astype(str),
                            'Muestra': adata_display.obs.get('sample', pd.NA).astype(str),
                            'N_Genes': adata_display.obs.get('n_genes_by_counts', pd.NA)
                        })
                        # ... (tu lógica para color_discrete_sequence_clusters si la tienes) ...
                        fig_umap_plotly_c = px.scatter(
                            df_umap_c_plotly, x='UMAP1', y='UMAP2', color='Cluster',
                            hover_data=['Muestra', 'N_Genes'],
                            title=f"UMAP Interactivo (Clusters Leiden, Res: {st.session_state.leiden_res})"
                        )
                        point_size_plotly_c = max(1, int(st.session_state.plot_point_size / 10))
                        fig_umap_plotly_c.update_traces(marker=dict(size=point_size_plotly_c, opacity=0.8))
                        st.plotly_chart(fig_umap_plotly_c, use_container_width=True)
                        # ... (botón de descarga HTML para el interactivo) ...
                    except Exception as e_plotly_c:
                        st.error(f"Error UMAP 2D interactivo (clústeres): {e_plotly_c}")
                        st.error(traceback.format_exc())

                # --- VERSIÓN ESTÁTICA (SCANPY/MATPLOTLIB) ---
                with st.expander("Ver UMAP 2D Estático por Clústeres (Scanpy)", expanded=False): # Contraído por defecto
                    try:
                        fig_umap_c_static, ax_c_static = plt.subplots(figsize=(7,6))
                        sc.pl.umap(
                            adata_display, 
                            color='leiden_clusters', 
                            legend_loc='on data', # <--- Muestra números/etiquetas en los clústeres
                            ax=ax_c_static, 
                            show=False, 
                            title=f"UMAP Estático (Clusters Leiden, Res: {st.session_state.leiden_res})", 
                            size=st.session_state.plot_point_size, 
                            palette=st.session_state.plot_palette if st.session_state.plot_palette != 'default' else None
                        )
                        st.pyplot(fig_umap_c_static)
                        st.download_button("UMAP Clústeres Estático (PNG)", fig_to_bytes(fig_umap_c_static), 
                                           "umap_clusters_static.png", "image/png", key="dl_umc_static_final")
                        plt.close(fig_umap_c_static)
                    except Exception as e_static_c:
                        st.error(f"Error UMAP 2D estático (clústeres): {e_static_c}")
                        st.error(traceback.format_exc())
            else: 
                st.warning("Clusters Leiden no encontrados para plot UMAP.")

            st.markdown("---") # Separador

            # --- UMAP 2D POR MUESTRA ---
            if 'sample' in adata_display.obs:
                st.subheader("UMAP 2D por Muestra")

                # --- VERSIÓN INTERACTIVA (PLOTLY) ---
                with st.expander("Ver UMAP 2D Interactivo por Muestra (Plotly)", expanded=True):
                    try:
                        # ... (tu código para UMAP 2D interactivo por Muestra con Plotly, como lo tenías) ...
                        umap_coords_s_plotly = adata_display.obsm['X_umap']
                        df_umap_s_plotly = pd.DataFrame({
                            'UMAP1': umap_coords_s_plotly[:, 0], 'UMAP2': umap_coords_s_plotly[:, 1],
                            'Muestra': adata_display.obs['sample'].astype(str),
                            'Cluster': adata_display.obs.get('leiden_clusters', "N/A").astype(str)
                        })
                        fig_umap_plotly_s = px.scatter(
                            df_umap_s_plotly, x='UMAP1', y='UMAP2', color='Muestra',
                            hover_data=['Cluster'], title="UMAP Interactivo por Muestra"
                        )
                        point_size_plotly_s = max(1, int(st.session_state.plot_point_size / 10))
                        fig_umap_plotly_s.update_traces(marker=dict(size=point_size_plotly_s, opacity=0.8))
                        st.plotly_chart(fig_umap_plotly_s, use_container_width=True)
                        # ... (botón de descarga HTML)
                    except Exception as e_plotly_s:
                        st.error(f"Error UMAP 2D interactivo (muestra): {e_plotly_s}")
                        st.error(traceback.format_exc())


                # --- VERSIÓN ESTÁTICA (SCANPY/MATPLOTLIB) ---
                with st.expander("Ver UMAP 2D Estático por Muestra (Scanpy)", expanded=False):
                    try:
                        fig_umap_s_static, ax_s_static = plt.subplots(figsize=(7,6))
                        sc.pl.umap(
                            adata_display, 
                            color='sample', 
                            ax=ax_s_static, 
                            show=False, 
                            title="UMAP Estático por Muestra", 
                            size=st.session_state.plot_point_size,
                            legend_loc='right margin', # O donde prefieras la leyenda de muestras
                            palette=st.session_state.plot_palette if st.session_state.plot_palette != 'default' else None
                        )
                        st.pyplot(fig_umap_s_static)
                        st.download_button("UMAP Muestra Estático (PNG)", fig_to_bytes(fig_umap_s_static), 
                                           "umap_sample_static.png", "image/png", key="dl_ums_static_final")
                        plt.close(fig_umap_s_static)
                    except Exception as e_static_s:
                        st.error(f"Error UMAP 2D estático (muestra): {e_static_s}")
                        st.error(traceback.format_exc())
            elif 'X_umap' in adata_display.obsm:
                 st.warning("Columna 'sample' no encontrada para el UMAP por Muestra.")

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
            st.subheader("UMAPs 2D por Muestra (Facetado, Coloreado por Clúster)")
            if 'sample' in adata_display.obs and 'leiden_clusters' in adata_display.obs and 'X_umap' in adata_display.obsm: # <--- CONDICIÓN 1
                try:
                    unique_samples_facet_disp = sorted(adata_display.obs['sample'].astype('category').cat.categories.tolist())
                    n_samples_facet_disp = len(unique_samples_facet_disp)
                    
                    if n_samples_facet_disp > 0: # <--- CONDICIÓN 2
                        cols_facet_disp = min(n_samples_facet_disp, 3)
                        rows_facet_disp = (n_samples_facet_disp + cols_facet_disp - 1) // cols_facet_disp
                        fig_facet_disp, axes_facet_disp = plt.subplots(rows_facet_disp, cols_facet_disp, 
                                                                    figsize=(cols_facet_disp * 5.5, rows_facet_disp * 5), squeeze=False)
                        axes_flat_facet_disp = axes_facet_disp.flatten()
                        idx_facet_disp = 0 # Inicializar por si el bucle no se ejecuta
                        
                        for idx_facet_disp, sample_val_facet_disp in enumerate(unique_samples_facet_disp):
                            if idx_facet_disp < len(axes_flat_facet_disp):
                                ax_curr_facet_disp = axes_flat_facet_disp[idx_facet_disp]
                                adata_subset_facet_disp = adata_display[adata_display.obs['sample'] == sample_val_facet_disp].copy() # <--- CREA SUBSET
                                
                                if not adata_subset_facet_disp.obs.empty and 'X_umap' in adata_subset_facet_disp.obsm: # <--- CONDICIÓN 3 (para el subset)
                                    sc.pl.umap(adata_subset_facet_disp, color='leiden_clusters', ax=ax_curr_facet_disp, show=False, 
                                               title=f"Muestra: {sample_val_facet_disp}", 
                                               legend_loc='on data' if idx_facet_disp == 0 and n_samples_facet_disp > 1 else None, 
                                               legend_fontsize=6, size=st.session_state.plot_point_size, 
                                               palette=st.session_state.plot_palette if st.session_state.plot_palette != 'default' else None)
                                elif not adata_subset_facet_disp.obs.empty: # Si hay células pero no X_umap en el subset
                                     ax_curr_facet_disp.text(0.5,0.5, f"M: {sample_val_facet_disp}\n(X_umap no disp.\nen subset)", ha='center',va='center', fontsize=8) # Mensaje más específico
                                else: # Si el subset está vacío
                                    ax_curr_facet_disp.text(0.5, 0.5, f"M: {sample_val_facet_disp}\n(Sin células\nen subset)", ha='center', va='center', fontsize=8)
                                ax_curr_facet_disp.set_xticks([]); ax_curr_facet_disp.set_yticks([]) # Limpiar ejes vacíos o con texto
                        
                        # Ocultar ejes no usados
                        for j_ax_empty_disp in range(idx_facet_disp + 1, len(axes_flat_facet_disp)): 
                            fig_facet_disp.delaxes(axes_flat_facet_disp[j_ax_empty_disp])
                        
                        plt.tight_layout()
                        st.pyplot(fig_facet_disp)
                        st.download_button("Descargar UMAPs Facetados (PNG)", fig_to_bytes(fig_facet_disp), "umaps_faceted.png", key="dl_umaps_facet_final_v2") # Nueva key
                        plt.close(fig_facet_disp)
                except Exception as e_facet_disp: 
                    st.error(f"Error UMAPs facetados: {e_facet_disp}")
                    st.error(traceback.format_exc()) # Mostrar traceback
            else: # Si falla la CONDICIÓN 1
                missing_keys_facet = []
                if 'sample' not in adata_display.obs: missing_keys_facet.append("'sample' en .obs")
                if 'leiden_clusters' not in adata_display.obs: missing_keys_facet.append("'leiden_clusters' en .obs")
                if 'X_umap' not in adata_display.obsm: missing_keys_facet.append("'X_umap' en .obsm")
                st.warning(f"No se pueden generar UMAPs facetados. Faltan datos necesarios: {', '.join(missing_keys_facet)}.")


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
                            print("DEBUG Heatmap: Dendrograma calculado. Ploteando CON dendrograma (método gcf).")
                            try:
                                # Crear una NUEVA figura para que Scanpy dibuje en ella
                                plt.figure(figsize=(10, max(6, len(genes_present_in_adata_for_heatmap) * 0.35)))
                                sc.pl.heatmap(
                                    adata_for_heatmap,
                                    genes_present_in_adata_for_heatmap,
                                    groupby='leiden_clusters',
                                    cmap=st.session_state.plot_palette if st.session_state.plot_palette != 'default' else None,
                                    standard_scale='var',
                                    dendrogram=True, # Debería encontrar la clave por defecto
                                    show=False
                                    # NO pasar 'ax' ni 'save'
                                )
                                fig_heatmap_actual = plt.gcf() # Obtener la figura actual
                                # plt.tight_layout() # Probar con y sin, por si causa problemas
                                st.pyplot(fig_heatmap_actual)
                                
                                bytes_img_heatmap = fig_to_bytes(fig_heatmap_actual)
                                st.download_button(
                                    "Descargar Heatmap Marcadores (PNG)",
                                    bytes_img_heatmap,
                                    "heatmap_markers.png", "image/png",
                                    key="dl_heatmap_markers_gcf_final_v3"
                                )
                                plt.close(fig_heatmap_actual)
                            except Exception as e_hm_gcf:
                                st.error(f"Error generando heatmap con método gcf: {e_hm_gcf}")
                                st.error(traceback.format_exc())


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
    

    with tab_gene_scoring_display:
        print("\n--- DEBUG: Inicio Pestaña Gene Scoring ---")
        print(f"adata_display.obs columns: {adata_display.obs.columns.tolist()}")
        print(f"adata_display.obsm keys: {list(adata_display.obsm.keys())}")
        #if score_to_visualize_selected in adata_display.obs: # Asumiendo que score_to_visualize_selected está definido
        #    print(f"Score '{score_to_visualize_selected}' dtype en Gene Scoring Tab: {adata_display.obs[score_to_visualize_selected].dtype}")
        st.subheader("Análisis de Scoring de Listas de Genes")
        st.write("""
        Esta sección permite calcular y visualizar un 'score' agregado para una o más listas de genes
        (firmas génicas) proporcionadas por el usuario. Es útil para investigar la actividad de programas
        génicos específicos o la presencia de tipos celulares definidos por conjuntos de marcadores.
        """)
        
        # Widgets para la entrada del usuario
        st.session_state.gene_score_user_lists_input = st.text_area(
            "Introduce tus listas de genes/firmas (una por línea):",
            value=st.session_state.gene_score_user_lists_input,
            height=200,
            key="gene_score_user_lists_ta_final_v2", # Key actualizada
            help="Formato: NombreFirma: GEN1,GEN2,...\nLos genes pueden estar separados por coma o espacio. Las líneas que empiezan con '#' son ignoradas."
        )
        st.session_state.gene_score_name_input = st.text_input( # Esta clave sí estaba en default_values
            "Nombre para la columna del Score (ej: LinfocitosT_Score):", 
            value=st.session_state.gene_score_name_input, # Usar la clave correcta de session_state
            key="gene_score_name_ti_final_v2" # Key del widget
        )

        if st.button("Calcular Scores de Firmas Génicas", key="calc_all_gene_scores_btn_v3"): # Key actualizada
            if not st.session_state.gene_score_user_lists_input.strip(): # Usar la clave correcta de session_state
                st.warning("Por favor, introduce al menos una lista de genes con su nombre.")
            else:
                parsed_gene_lists_gs = {} # Renombrar para evitar conflicto
                found_any_valid_genes_gs = False

                for line_gs in st.session_state.gene_score_user_lists_input.splitlines(): # Usar la clave correcta
                    line_gs = line_gs.strip()
                    if not line_gs or line_gs.startswith('#'): continue
                    
                    parts_gs = line_gs.split(':', 1)
                    if len(parts_gs) != 2:
                        st.error(f"Línea mal formateada: '{line_gs}'.\nFormato esperado: 'NombreDeLaFirma: GENE1, GENE2,...'\nAsegúrate de incluir un nombre para la firma seguido de ':' y luego la lista de genes.")
                        continue
                    
                    score_name_raw_gs = parts_gs[0].strip()
                    score_name_col_gs = "".join(c if c.isalnum() or c == '_' else '_' for c in score_name_raw_gs).strip('_')
                    if not score_name_col_gs: 
                        score_name_col_gs = f"score_{pd.Timestamp.now().strftime('%H%M%S%f')}" # Más único
                        st.warning(f"Nombre de firma '{score_name_raw_gs}' inválido, usando '{score_name_col_gs}'.")

                    genes_str_gs = parts_gs[1].strip()
                    genes_in_list_raw_gs = [g.strip() for g in genes_str_gs.replace(',', ' ').split() if g.strip()]

                    if not genes_in_list_raw_gs:
                        st.warning(f"La firma '{score_name_raw_gs}' no contiene genes.")
                        continue

                    current_valid_genes_gs = []
                    current_not_found_gs = []
                    # valid_genes_lower_map_display debe estar definido si estamos dentro del if principal
                    if not valid_genes_lower_map_display:
                         st.error("Mapa de genes no disponible para validación.") # No debería ocurrir aquí
                    else:
                        for g_raw_gs in genes_in_list_raw_gs:
                            g_lower_gs = g_raw_gs.lower()
                            if g_lower_gs in valid_genes_lower_map_display:
                                current_valid_genes_gs.append(valid_genes_lower_map_display[g_lower_gs])
                            else:
                                current_not_found_gs.append(g_raw_gs)
                    
                    if current_not_found_gs:
                        st.warning(f"Para firma '{score_name_raw_gs}': Genes no encontrados -> {', '.join(current_not_found_gs)}")
                    
                    if current_valid_genes_gs:
                        parsed_gene_lists_gs[score_name_col_gs] = current_valid_genes_gs
                        found_any_valid_genes_gs = True
                    else:
                        st.error(f"Ningún gen válido para firma '{score_name_raw_gs}'.")
                
                if found_any_valid_genes_gs:
                    with st.spinner("Calculando scores..."):
                        adata_for_scoring_gs = st.session_state.adata_processed # Trabajar sobre el procesado
                        
                        for score_col_loop, gene_l_loop in parsed_gene_lists_gs.items():
                            try:
                                sc.tl.score_genes(
                                    adata_for_scoring_gs, 
                                    gene_list=gene_l_loop, 
                                    score_name=score_col_loop, 
                                    random_state=0, use_raw=False
                                )
                                st.success(f"Score '{score_col_loop}' calculado y añadido a `adata.obs`.")
                                if score_col_loop not in st.session_state.gene_scores_calculated:
                                    st.session_state.gene_scores_calculated[score_col_loop] = "Calculado"
                                # Actualizar explícitamente el AnnData en session_state
                                st.session_state.adata_processed = adata_for_scoring_gs 
                            except Exception as e_score_calc_gs:
                                st.error(f"Error calculando score para '{score_col_loop}': {e_score_calc_gs}")
                                st.session_state.gene_scores_calculated[score_col_loop] = f"Error: {e_score_calc_gs}"
                        st.rerun() 
                else:
                    st.error("No se encontraron listas de genes válidas para calcular scores.")

        # --- Visualización de los scores calculados ---
        st.markdown("---")
        st.subheader("Visualizar Scores de Firmas Calculadas")
        
        calculated_score_options = [
            sc_name for sc_name, status in st.session_state.get("gene_scores_calculated", {}).items() 
            if status == "Calculado" and sc_name in adata_display.obs.columns
        ]

        if not calculated_score_options:
            st.info("No hay scores calculados disponibles para visualizar. Por favor, calcula un score primero.")
        else:
            score_name_from_input_gs = st.session_state.gene_score_name_input.strip() # Nombre del último text input
            default_selection_idx_gs = 0
            if score_name_from_input_gs in calculated_score_options:
                default_selection_idx_gs = calculated_score_options.index(score_name_from_input_gs)

            score_to_visualize_selected = st.selectbox(
                "Selecciona un score calculado para visualizar:", 
                options=calculated_score_options,
                index=default_selection_idx_gs,
                key="select_score_to_viz_final_v3" # Nueva key
            )

            if score_to_visualize_selected:
                st.markdown(f"#### Visualización del Score: **{score_to_visualize_selected}**")
                
                # --- DEBUG ANTES DE PLOTS DE SCORE ---
                print(f"\n--- DEBUG: Visualizando Score: {score_to_visualize_selected} ---")
                if score_to_visualize_selected not in adata_display.obs.columns:
                    st.error(f"Error interno: El score seleccionado '{score_to_visualize_selected}' no se encuentra en adata_display.obs.")
                    print(f"ERROR: Score '{score_to_visualize_selected}' no en adata_display.obs. Columnas: {adata_display.obs.columns.tolist()}")
                    # No continuar con los plots si el score no existe
                else:
                    print(f"Valores del score '{score_to_visualize_selected}' (primeros 5): {adata_display.obs[score_to_visualize_selected].head().tolist()}")
                    print(f"Tipo de datos del score: {adata_display.obs[score_to_visualize_selected].dtype}")
                    print(f"Hay NaNs en el score?: {adata_display.obs[score_to_visualize_selected].isnull().any()}")
                    print(f"Hay Infs en el score?: {np.isinf(adata_display.obs[score_to_visualize_selected]).any()}")
                
                    col_score_viz1_disp_gs, col_score_viz2_disp_gs = st.columns(2) # Nombres de variables únicos
                    
                    with col_score_viz1_disp_gs:
                        if 'X_umap' in adata_display.obsm:
                            st.markdown("##### Score en UMAP 2D")
                            print(f"DEBUG UMAP 2D Score: Intentando plotear score '{score_to_visualize_selected}'")
                            fig_s_u2d_gs_plot = None # Inicializar por si falla la creación
                            try:
                                # Crear figura y eje
                                fig_s_u2d_gs_plot, ax_s_u2d_gs_plot = plt.subplots(figsize=(7,6)) # Un figsize un poco más grande

                                # Llamada simplificada a sc.pl.umap
                                sc.pl.umap(
                                    adata_display, 
                                    color=score_to_visualize_selected, 
                                    ax=ax_s_u2d_gs_plot, 
                                    show=False,
                                    legend_loc=None, # Quitar leyenda si no es necesaria para un score continuo
                                    cmap='viridis' # Forzar un cmap conocido
                                    # Quitar size y palette temporalmente
                                )
                                ax_s_u2d_gs_plot.set_title(f"Score: {score_to_visualize_selected}", fontsize=12)
                                
                                # Intentar st.pyplot ANTES de tight_layout
                                st.pyplot(fig_s_u2d_gs_plot) 
                                print(f"DEBUG UMAP 2D Score: st.pyplot para '{score_to_visualize_selected}' ejecutado.")

                                # plt.tight_layout() # Comentar temporalmente si da warnings
                                
                            except Exception as e_sup2d_detailed:
                                st.error(f"Error DETALLADO UMAP 2D para score '{score_to_visualize_selected}': {e_sup2d_detailed}")
                                st.error(traceback.format_exc()) # Mostrar el traceback completo
                            finally:
                                if fig_s_u2d_gs_plot is not None: # Solo cerrar si la figura fue creada
                                    plt.close(fig_s_u2d_gs_plot)
                                    print(f"DEBUG UMAP 2D Score: Figura para '{score_to_visualize_selected}' cerrada.")
                        else:
                            st.warning("Coordenadas UMAP 2D no disponibles para plotear score.")
                    
                    with col_score_viz2_disp_gs:
                        if 'leiden_clusters' in adata_display.obs:
                            st.markdown("##### Score en Violín por Clúster")
                            print(f"DEBUG Violín Score: Intentando plotear score '{score_to_visualize_selected}'")
                            fig_s_vln_gs_plot = None # Inicializar
                            try:
                                n_clusters_vln_score = adata_display.obs['leiden_clusters'].nunique()
                                fig_s_vln_gs_plot, ax_s_vln_gs_plot = plt.subplots(figsize=(max(6, n_clusters_vln_score * 0.8), 5))

                                sc.pl.violin(
                                    adata_display, 
                                    keys=score_to_visualize_selected, 
                                    groupby='leiden_clusters', 
                                    rotation=45, 
                                    ax=ax_s_vln_gs_plot, 
                                    show=False,
                                    use_raw=False,
                                    cut=0
                                    # Quitar title de aquí
                                )
                                ax_s_vln_gs_plot.set_title(f"Score: {score_to_visualize_selected} por Clúster", fontsize=10)
                                
                                st.pyplot(fig_s_vln_gs_plot)
                                print(f"DEBUG Violín Score: st.pyplot para '{score_to_visualize_selected}' ejecutado.")
                                
                                # plt.tight_layout() # Comentar temporalmente

                            except Exception as e_svln_detailed:
                                st.error(f"Error DETALLADO Violín para score '{score_to_visualize_selected}': {e_svln_detailed}")
                                st.error(traceback.format_exc())
                            finally:
                                if fig_s_vln_gs_plot is not None:
                                    plt.close(fig_s_vln_gs_plot)
                                    print(f"DEBUG Violín Score: Figura para '{score_to_visualize_selected}' cerrada.")
                        else:
                            st.warning("Clusters Leiden no disponibles para plotear score en violín.")

                    if st.session_state.calc_umap_3d:
                        if 'X_umap_3d' in adata_display.obsm:
                            st.markdown(f"##### Score '{score_to_visualize_selected}' en UMAP 3D")
                            try:
                                umap_3d_coords_score_plot = adata_display.obsm['X_umap_3d']
                                df_umap3d_gene_score_plot = pd.DataFrame({ # Nombres de variables únicos
                                    'UMAP1': umap_3d_coords_score_plot[:, 0],
                                    'UMAP2': umap_3d_coords_score_plot[:, 1],
                                    'UMAP3': umap_3d_coords_score_plot[:, 2],
                                    score_to_visualize_selected: adata_display.obs[score_to_visualize_selected],
                                    'Cluster': adata_display.obs.get('leiden_clusters', pd.NA).astype(str),
                                    'Muestra': adata_display.obs.get('sample', pd.NA).astype(str)
                                })
                                n_clusters_3d_score_plot = adata_display.obs.get('leiden_clusters', pd.Series(dtype=str)).nunique()
                                color_seq_3d_score = px.colors.qualitative.Plotly if st.session_state.plot_palette == "default" or n_clusters_3d_score_plot > len(px.colors.qualitative.Plotly) else getattr(px.colors.qualitative, st.session_state.plot_palette, px.colors.qualitative.Plotly)
                                
                                fig_3d_gene_score_plot = px.scatter_3d( # Nombres de variables únicos
                                    df_umap3d_gene_score_plot, x='UMAP1', y='UMAP2', z='UMAP3',
                                    color=score_to_visualize_selected, color_continuous_scale='viridis', # Usar viridis para score continuo
                                    hover_data=['Muestra', 'Cluster', score_to_visualize_selected],
                                    title=f"UMAP 3D - Score: {score_to_visualize_selected}"
                                )
                                marker_size_3d_score_plot = max(1, int(st.session_state.plot_point_size / 15))
                                fig_3d_gene_score_plot.update_traces(marker=dict(size=marker_size_3d_score_plot))
                                st.plotly_chart(fig_3d_gene_score_plot, use_container_width=True)
                            except Exception as e_plot3d_score_tab_plot: # Nombres de variables únicos
                                st.error(f"Error UMAP 3D para score '{score_to_visualize_selected}': {e_plot3d_score_tab_plot}")
                                st.error(traceback.format_exc())
                        elif score_to_visualize_selected: # Si calc_umap_3d es True pero X_umap_3d no está
                            st.info(f"UMAP 3D fue seleccionado pero sus coordenadas no están disponibles para visualizar el score '{score_to_visualize_selected}'.")
            # No se necesita un 'else' aquí porque el selectbox siempre tendrá un valor si calculated_score_options no está vacío.
            
            elif st.session_state.gene_score_user_lists_input.strip() and score_name_from_input_gs: # score_name_from_input_gs es el del text_input
                # Este elif se activa si el usuario ha escrito algo en los inputs pero NO hay scores calculados válidos.
                # O si el score que está en el input text no se pudo calcular.
                if score_name_from_input_gs not in calculated_score_options: # Y no está entre los calculados con éxito
                    st.info(f"El score '{score_name_from_input_gs}' no se encuentra en los datos o no se pudo calcular. Verifica los mensajes de error/warning arriba o calcúlalo.")
 
    with tab_qc_display:
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

    with tab_dea_display:
        st.subheader("Resultados del Análisis de Expresión Diferencial")
        if st.session_state.dea_results_df is not None and not st.session_state.dea_results_df.empty:
            st.markdown(f"**Comparación Actual:** `{st.session_state.dea_comparison_str}`")
            
            # --- Configuración de Columnas para st.data_editor ---
            column_config_dea = {
                "P-Value": st.column_config.NumberColumn(
                    "P-Valor",
                    format="%.2e", # Notación científica
                    help="P-valor crudo."
                ),
                "P-Value Adj": st.column_config.NumberColumn(
                    "P-Valor Ajustado",
                    format="%.2e", # Notación científica
                    help="P-valor ajustado por múltiples comparaciones."
                ),
                "Scores": st.column_config.NumberColumn( # 'Scores' en plural como en tu DataFrame
                    "Score",
                    format="%.3f"
                ),
                "Log2FC": st.column_config.NumberColumn(
                    "Log2 Fold Change",
                    format="%.3f"
                ),
                "Gene": st.column_config.TextColumn("Gen")
            }

            st.data_editor(
                st.session_state.dea_results_df.head(st.session_state.dea_n_genes_display), # Muestra N genes
                height=400, 
                use_container_width=True, 
                num_rows="dynamic", # Para ver más de los N genes si se desea
                column_config=column_config_dea,
                key="dea_results_data_editor_final_v2" # Key nueva
            )
        elif st.session_state.analysis_done:
            st.info("No hay resultados de DEA. Ejecuta el análisis desde la sidebar si es necesario.")


    with tab_gene_explorer_display:
        st.subheader("Visualización de Expresión para Genes Específicos")
        if not genes_to_visualize_list: 
            st.info("Ingresa nombres de genes válidos para visualizarlos.")
        else:
            st.write(f"Mostrando expresión para: **{', '.join(genes_to_visualize_list)}**")
            # UMAPs por Expresión Génica
            if 'X_umap' not in adata_display.obsm:
             st.warning("Plots UMAP no disponibles (UMAP no calculado o falló).")
            else:
                st.markdown("#### UMAPs coloreados por Expresión Génica")
                if genes_to_visualize_list: # Si hay genes en la lista
                    selected_gene_for_interactive_umap = st.selectbox(
                        "Selecciona un gen para UMAP Interactivo:", 
                        options=genes_to_visualize_list,
                        key="select_gene_interactive_umap"
                    )
                    if selected_gene_for_interactive_umap:
                        try:
                            umap_coords_ge = adata_display.obsm['X_umap']
                            df_umap_ge = pd.DataFrame({
                                'UMAP1': umap_coords_ge[:, 0],
                                'UMAP2': umap_coords_ge[:, 1],
                                'Expresión': adata_display[:, selected_gene_for_interactive_umap].X.toarray().flatten(), # Asegurar que es denso y 1D
                                'Cluster': adata_display.obs.get('leiden_clusters', pd.NA).astype(str)
                            })
                            fig_umap_plotly_ge = px.scatter(
                                df_umap_ge, x='UMAP1', y='UMAP2', color='Expresión',
                                color_continuous_scale='viridis', # Para expresión continua
                                hover_data=['Cluster'],
                                title=f"UMAP 2D Interactivo: Expresión de {selected_gene_for_interactive_umap}"
                            )
                            point_size_plotly_ge = max(1, int(st.session_state.plot_point_size / 10))
                            fig_umap_plotly_ge.update_traces(marker=dict(size=point_size_plotly_ge, opacity=0.8))
                            st.plotly_chart(fig_umap_plotly_ge, use_container_width=True)
                            # ... (botón de descarga HTML) ...
                        except Exception as e_ge_plotly_umap:
                            st.error(f"Error UMAP interactivo para gen {selected_gene_for_interactive_umap}: {e_ge_plotly_umap}")

            st.markdown("---")
            st.markdown("##### UMAPs Estáticos (Múltiples Genes)")             
            if 'X_umap' not in adata_display.obsm:
                st.warning("Plots UMAP no disponibles (UMAP no calculado o falló).")
            else:
                st.markdown("#### UMAPs coloreados por Expresión Génica")
                n_genes_umap_plot_exp = len(genes_to_visualize_list) # Usar la variable correcta
                cols_genes_umap_exp = min(n_genes_umap_plot_exp, 3)
                rows_genes_umap_exp = (n_genes_umap_plot_exp + cols_genes_umap_exp - 1) // cols_genes_umap_exp
                
                if n_genes_umap_plot_exp > 0:
                    fig_ge_umaps_exp, axes_ge_umaps_exp = plt.subplots(
                        rows_genes_umap_exp, 
                        cols_genes_umap_exp, 
                        figsize=(cols_genes_umap_exp * 5, rows_genes_umap_exp * 4.5), 
                        squeeze=False
                    )
                    axes_flat_ge_umaps_exp = axes_ge_umaps_exp.flatten()
                    idx_ge_plot_exp = 0 
                    for idx_ge_plot_exp, gene_name_plot_exp in enumerate(genes_to_visualize_list): # Usar la variable correcta
                        if idx_ge_plot_exp < len(axes_flat_ge_umaps_exp):
                            ax_ge_curr_exp = axes_flat_ge_umaps_exp[idx_ge_plot_exp]
                            try:
                                sc.pl.umap(adata_display, color=gene_name_plot_exp, ax=ax_ge_curr_exp, 
                                           show=False, title=gene_name_plot_exp, cmap='viridis', 
                                           use_raw=False, size=st.session_state.plot_point_size)
                            except Exception as e_ge_umap_plot_exp: 
                                ax_ge_curr_exp.text(0.5, 0.5, f"Error plot\n{gene_name_plot_exp}", ha='center', va='center', color='red')
                                ax_ge_curr_exp.set_xticks([]); ax_ge_curr_exp.set_yticks([])
                                print(f"Error ploteando UMAP para gen {gene_name_plot_exp}: {e_ge_umap_plot_exp}") 
                    
                    for j_ge_empty_ax_exp in range(idx_ge_plot_exp + 1, len(axes_flat_ge_umaps_exp)): 
                        fig_ge_umaps_exp.delaxes(axes_flat_ge_umaps_exp[j_ge_empty_ax_exp])
                    
                    plt.tight_layout()
                    st.pyplot(fig_ge_umaps_exp)
                    st.download_button("Descargar UMAPs de Genes (PNG)", fig_to_bytes(fig_ge_umaps_exp), 
                                       "gene_explorer_umaps.png", "image/png", key="dl_ge_umaps_png_final") # Key actualizada
                    plt.close(fig_ge_umaps_exp)
            
            # --- VIOLINES POR CLÚSTER (CORREGIDO) ---
            if 'leiden_clusters' in adata_display.obs and genes_to_visualize_list: # Usar la variable correcta
                st.markdown("#### Diagramas de Violín por Clúster de Leiden")
                
                n_genes_vln_cl_exp = len(genes_to_visualize_list) # Usar la variable correcta
                cols_vln_cl_exp = min(2, n_genes_vln_cl_exp) 
                rows_vln_cl_exp = (n_genes_vln_cl_exp + cols_vln_cl_exp - 1) // cols_vln_cl_exp

                if n_genes_vln_cl_exp > 0:
                    fig_violins_cl_exp, axes_violins_cl_exp = plt.subplots(
                        rows_vln_cl_exp, cols_vln_cl_exp, 
                        figsize=(cols_vln_cl_exp * max(7, adata_display.obs['leiden_clusters'].nunique() * 0.6), rows_vln_cl_exp * 5), 
                        squeeze=False 
                    )
                    axes_flat_cl_exp = axes_violins_cl_exp.flatten()
                    
                    idx_plot_vln_cl_actual = 0 
                    for idx_plot_vln_cl_actual, gene_name_vln_cl_exp in enumerate(genes_to_visualize_list): # Usar la variable correcta
                        if idx_plot_vln_cl_actual < len(axes_flat_cl_exp): 
                            ax_curr_vln_cl_exp = axes_flat_cl_exp[idx_plot_vln_cl_actual]
                            try:
                                sc.pl.violin(
                                    adata_display, 
                                    keys=gene_name_vln_cl_exp, # <--- Variable del bucle
                                    groupby='leiden_clusters', 
                                    rotation=45,
                                    ax=ax_curr_vln_cl_exp, 
                                    show=False, 
                                    use_raw=False, 
                                    cut=0
                                    # No title aquí
                                )
                                ax_curr_vln_cl_exp.set_title(gene_name_vln_cl_exp) # <--- ¡ERROR AQUÍ! Debe ser gene_name_vln_cl_exp
                            except Exception as e_vln_gene_cl_exp_loop: 
                                ax_curr_vln_cl_exp.text(0.5,0.5, f"Error plot\n{gene_name_vln_cl_exp}", ha='center', va='center', color='red')
                                ax_curr_vln_cl_exp.set_xticks([]); ax_curr_vln_cl_exp.set_yticks([])
                                print(f"Error al generar violín (clúster) para gen '{gene_name_vln_cl_exp}': {e_vln_gene_cl_exp_loop}")
                        else: break
                    
                    for j_empty_ax_vln_cl_exp in range(idx_plot_vln_cl_actual + 1, len(axes_flat_cl_exp)):
                        fig_violins_cl_exp.delaxes(axes_flat_cl_exp[j_empty_ax_vln_cl_exp])

                    plt.tight_layout()
                    st.pyplot(fig_violins_cl_exp)
                    st.download_button("Descargar Violines por Clúster (PNG)", fig_to_bytes(fig_violins_cl_exp), 
                                       "ge_violins_cluster.png", key="dl_ge_violins_cluster_final_v2") 
                    plt.close(fig_violins_cl_exp)
            
            # --- VIOLINES POR CONDICIÓN (CORREGIDO) ---
            if 'condition_temp_dea' in adata_display.obs and adata_display.obs['condition_temp_dea'].nunique() > 1 and genes_to_visualize_list: # Usar la variable correcta
                st.markdown("#### Diagramas de Violín por Condición (definida en DEA)")
                
                n_genes_vln_cond_exp = len(genes_to_visualize_list) # Usar la variable correcta
                cols_vln_cond_exp = min(2, n_genes_vln_cond_exp)
                rows_vln_cond_exp = (n_genes_vln_cond_exp + cols_vln_cond_exp - 1) // cols_vln_cond_exp

                if n_genes_vln_cond_exp > 0:
                    fig_violins_cond_exp, axes_violins_cond_exp = plt.subplots(
                        rows_vln_cond_exp, cols_vln_cond_exp, 
                        figsize=(cols_vln_cond_exp * max(7, adata_display.obs['condition_temp_dea'].nunique() * 1.0), rows_vln_cond_exp * 5),
                        squeeze=False
                    )
                    axes_flat_cond_exp = axes_violins_cond_exp.flatten()
                    
                    idx_plot_vln_cond_actual = 0
                    for idx_plot_vln_cond_actual, gene_name_vln_cond_exp in enumerate(genes_to_visualize_list): # Usar la variable correcta
                        if idx_plot_vln_cond_actual < len(axes_flat_cond_exp):
                            ax_curr_vln_cond_exp = axes_flat_cond_exp[idx_plot_vln_cond_actual]
                            try:
                                sc.pl.violin(adata_display, keys=gene_name_vln_cond_exp, groupby='condition_temp_dea', rotation=45,
                                             ax=ax_curr_vln_cond_exp, show=False, use_raw=False, cut=0, title=gene_name_vln_cond_exp)
                            except Exception as e_vln_gene_cond_exp_loop:
                                ax_curr_vln_cond_exp.text(0.5,0.5, f"Error plot\n{gene_name_vln_cond_exp}", ha='center', va='center', color='red')
                                ax_curr_vln_cond_exp.set_xticks([]); ax_curr_vln_cond_exp.set_yticks([])
                                print(f"Error al generar violín (condición) para gen '{gene_name_vln_cond_exp}': {e_vln_gene_cond_exp_loop}")
                        else: break
                    
                    for j_empty_ax_vln_cond_exp in range(idx_plot_vln_cond_actual + 1, len(axes_flat_cond_exp)):
                        fig_violins_cond_exp.delaxes(axes_flat_cond_exp[j_empty_ax_vln_cond_exp])

                    plt.tight_layout()
                    st.pyplot(fig_violins_cond_exp)
                    st.download_button("Descargar Violines por Condición (PNG)", fig_to_bytes(fig_violins_cond_exp), 
                                       "ge_violins_condition.png", key="dl_ge_violins_condition_final_v2")
                    plt.close(fig_violins_cond_exp)

            # --- DOT PLOT DEL EXPLORADOR DE GENES ---
            if len(genes_to_visualize_list) > 0 and 'leiden_clusters' in adata_display.obs: # Usar la variable correcta
                st.markdown("#### Dot Plot de Genes Seleccionados por Clúster")
                try:
                    n_clusters_dot_ge_exp = adata_display.obs['leiden_clusters'].nunique()
                    fig_ge_dotplot_exp, ax_ge_dotplot_exp = plt.subplots(figsize=(max(8, len(genes_to_visualize_list) * 0.7), max(5, n_clusters_dot_ge_exp * 0.5))) # Usar la variable correcta
                    sc.pl.dotplot(adata_display, genes_to_visualize_list, groupby='leiden_clusters', ax=ax_ge_dotplot_exp, show=False, standard_scale='var', use_raw=False) # Usar la variable correcta
                    plt.xticks(rotation=90); plt.tight_layout(); st.pyplot(fig_ge_dotplot_exp)
                    st.download_button("Descargar Dot Plot de Genes (PNG)", fig_to_bytes(fig_ge_dotplot_exp), "ge_dotplot.png", key="dl_ge_dotplot_final") # Key actualizada
                    plt.close(fig_ge_dotplot_exp)
                except Exception as e_ge_dot_exp: 
                    st.error(f"Error dot plot genes seleccionados: {e_ge_dot_exp}")
                    st.error(traceback.format_exc()) # Añadir traceback para depuración
                    st.info("No se pudo generar el dot plot. Asegúrate de que los genes están en el dataset y los clusters están definidos.")
            else:
                st.info("No se encontraron genes seleccionados o clusters definidos para el dot plot.")


    with tab_info_display: # O como hayas llamado a la variable de tu pestaña de Info
        st.subheader("Información del Dataset Procesado y Diagnóstico PCA") # Modifica el subheader si quieres
        
        # --- INFORMACIÓN BÁSICA DEL DATASET (COMO LA TENÍAS) ---
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
        
        # --- INICIO: PLOT DE VARIANZA PCA (ELBOW PLOT) ---
        if st.session_state.get("show_pca_variance", True): # Usar .get para default si la clave no existiera
            st.markdown("---") 
            st.markdown("#### Diagnóstico de Componentes Principales (PCA)")
            
            # Usar el adata_hvg_subset donde se calculó PCA
            adata_for_pca_diagnosis = st.session_state.adata_hvg_subset 
            
            if adata_for_pca_diagnosis is not None and \
               'pca' in adata_for_pca_diagnosis.uns and \
               'variance_ratio' in adata_for_pca_diagnosis.uns['pca']:
                
                variance_ratio_data_pca = adata_for_pca_diagnosis.uns['pca']['variance_ratio']
                n_pcs_calculated_pca = len(variance_ratio_data_pca)
                
                # Obtener el número de PCs que el pipeline realmente usó
                pcs_actually_used_pca = st.session_state.get("n_pcs_actually_used_in_pipeline", st.session_state.n_pcs)

                col_elbow_pca, col_cumvar_pca = st.columns(2)

                with col_elbow_pca:
                    st.markdown("##### Varianza Explicada por cada PC")
                    fig_elbow_plot_pca, ax_elbow_plot_pca = plt.subplots(figsize=(7,4))
                    ax_elbow_plot_pca.plot(range(1, n_pcs_calculated_pca + 1), variance_ratio_data_pca, marker='o', linestyle='-', color='dodgerblue')
                    ax_elbow_plot_pca.set_xlabel("Componente Principal")
                    ax_elbow_plot_pca.set_ylabel("Proporción de Varianza Explicada")
                    ax_elbow_plot_pca.set_title("Elbow Plot para Selección de PCs")
                    ax_elbow_plot_pca.grid(True, linestyle=':', alpha=0.7)
                    ax_elbow_plot_pca.axvline(x=pcs_actually_used_pca, color='red', linestyle='--', 
                                              label=f'PCs Usados: {pcs_actually_used_pca}')
                    ax_elbow_plot_pca.legend()
                    st.pyplot(fig_elbow_plot_pca)
                    plt.close(fig_elbow_plot_pca)

                with col_cumvar_pca:
                    st.markdown("##### Varianza Acumulada Explicada")
                    cumulative_variance_pca = np.cumsum(variance_ratio_data_pca)
                    fig_cumvar_plot_pca, ax_cumvar_plot_pca = plt.subplots(figsize=(7,4))
                    ax_cumvar_plot_pca.plot(range(1, n_pcs_calculated_pca + 1), cumulative_variance_pca, marker='.', linestyle='-', color='orangered')
                    ax_cumvar_plot_pca.set_xlabel("Número de Componentes Principales")
                    ax_cumvar_plot_pca.set_ylabel("Varianza Acumulada Explicada")
                    ax_cumvar_plot_pca.set_title("Varianza Acumulada de PCA")
                    ax_cumvar_plot_pca.grid(True, linestyle=':', alpha=0.7)
                    ax_cumvar_plot_pca.axhline(y=0.9, color='grey', linestyle=':', label='90% Varianza')
                    ax_cumvar_plot_pca.axhline(y=0.8, color='lightgrey', linestyle=':', label='80% Varianza')
                    ax_cumvar_plot_pca.axvline(x=pcs_actually_used_pca, color='red', linestyle='--', 
                                               label=f'PCs Usados: {pcs_actually_used_pca}')
                    ax_cumvar_plot_pca.legend()
                    st.pyplot(fig_cumvar_plot_pca)
                    plt.close(fig_cumvar_plot_pca)

                # Mostrar tabla de varianza (opcional)
                # df_pca_variance_table_info = pd.DataFrame({
                #     'PC': range(1, n_pcs_calculated_pca + 1),
                #     'Varianza Individual': variance_ratio_data_pca,
                #     'Varianza Acumulada': cumulative_variance_pca
                # })
                # st.markdown("##### Tabla de Varianza por PC (primeros PCs)")
                # st.dataframe(df_pca_variance_table_info.head(max(20, pcs_actually_used_pca + 5)))
            else:
                st.info("Datos de varianza PCA no disponibles. Ejecuta el pipeline principal y asegúrate de que el PCA se calcule sobre 'adata_hvg_subset'.")
        elif not st.session_state.get("show_pca_variance", True) and st.session_state.analysis_done : # Si el análisis se hizo pero el usuario desactivó el plot
             st.info("La visualización de la varianza PCA está desactivada en la configuración de la sidebar ('Personalización de Plots').")
        # --- FIN: PLOT DE VARIANZA PCA ---

        st.markdown("---") # Separador
        st.write("Primeras 5 filas de Metadatos de Células (`.obs`):")
        st.dataframe(adata_display.obs.head())
        st.write("Primeras 5 filas de Metadatos de Genes (`.var`):")
        st.dataframe(adata_display.var.head())
 
        st.markdown("---") # Separador
        st.write("Si necesitas más información, revisa el objeto `adata` completo en la consola de Streamlit.")
        st.info("Recuerda que puedes descargar los resultados de cada paso del pipeline desde la sidebar.")

else: # Si el análisis no se ha completado
    # Recalcular all_files_provided aquí para este scope
    all_files_provided_main_scope = True
    if st.session_state.num_samples > 0:
        for i in range(st.session_state.num_samples):
            if not (st.session_state.sample_files.get(f"matrix_file_{i}") and \
                    st.session_state.sample_files.get(f"features_file_{i}") and \
                    st.session_state.sample_files.get(f"barcodes_file_{i}")):
                all_files_provided_main_scope = False; break
    else: all_files_provided_main_scope = False

    if st.session_state.adata_raw is None and st.session_state.num_samples > 0 and not all_files_provided_main_scope:
        st.info("Bienvenido. Sube todos los archivos y haz clic en 'Cargar y Concatenar Datos'.")
    elif st.session_state.num_samples < 1 : 
        st.info("Bienvenido. Ajusta el 'Número de muestras' (mínimo 1) para comenzar.")
    elif st.session_state.adata_raw is None and st.session_state.num_samples > 0 and all_files_provided_main_scope:
        st.info("Archivos listos. Haz clic en 'Cargar y Concatenar Datos' en la sidebar.")
    elif st.session_state.adata_raw is not None and not st.session_state.analysis_done:
        st.info("Datos cargados. Haz clic en 'Ejecutar Pipeline Principal' en la sidebar.")
    else: 
        st.info("Bienvenido al Analizador Interactivo de scRNA-seq. Configura tus muestras y parámetros en la sidebar para comenzar.")

# Notas finales en la Sidebar
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


# Nota sobre dependencias y versiones
st.sidebar.markdown("---")
st.sidebar.info("Analizador scRNA-seq v1.0. Basado en Scanpy.")
st.sidebar.markdown("Si experimentas errores con UMAP (ej: `ValueError: high is out of bounds`), considera usar un entorno con `numpy<2.0` o Python 3.10/3.11.")
st.sidebar.markdown("Se recomienda crear un entorno virtual con versiones compatibles de las bibliotecas (ej: `numpy<2.0` si se experimentan errores con UMAP).")
st.sidebar.markdown("Si tienes problemas, consulta la [documentación de Scanpy](https://scanpy.readthedocs.io/en/stable/) o el [repositorio de GitHub]).")
st.sidebar.markdown("Para más información, visita el [repositorio de GitHub]).")
st.sidebar.markdown("**Desarrollado por:** Pedro Botías - pbotias@ucm.es - https://github.com/pbotiast/scRNASeq")
st.sidebar.markdown("**Licencia:** Licencia MIT - https://opensource.org/licenses/MIT")
st.sidebar.markdown("**Fecha:** [05/05/2025] - [20/05/2025]")  
st.sidebar.markdown("**Versión:** 1.1")
st.sidebar.markdown("**Última Actualización:** 2025-05-20")
st.sidebar.markdown("**Notas:** Esta aplicación es un prototipo y puede contener errores. Usa bajo tu propio riesgo.")
st.sidebar.markdown("**Disclaimer:** Esta aplicación es un prototipo y puede contener errores. Usa bajo tu propio riesgo.")








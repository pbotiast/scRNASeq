# -*- coding: utf-8 -*-
# Copyright (c) 2025 Pedro Botías
# Licenciado bajo la Licencia MIT; consulta LICENSE.txt para más detalles.
"""
Streamlit app para análisis interactivo de datos de Single-Cell RNA-seq con múltiples muestras.
Permite cargar datos 10x, realizar QC, normalización, PCA, UMAP, clustering, DEA y explorar genes.
Versión mejorada con optimización de memoria, manejo de errores y UX mejorada.
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
import difflib
import gzip
from scipy import sparse
# import warnings # No warnings.catch_warnings activos por ahora

# Configuración de la página
st.set_page_config(layout="wide", initial_sidebar_state="expanded")
st.title("Analizador Interactivo de Single-Cell RNA-seq")

# --- Funciones Auxiliares ---
def validate_10x_files(matrix_file, features_file, barcodes_file):
    """Valida que los archivos 10x tengan el formato correcto."""
    files_to_check = {
        "Matrix": (matrix_file, "%%MatrixMarket"),
        "Features": (features_file, None), # Para TSV, solo intentamos leer
        "Barcodes": (barcodes_file, None)  # Para TSV, solo intentamos leer
    }
    try:
        for name, (file_obj, expected_start) in files_to_check.items():
            if file_obj is None: # Si el archivo no fue subido
                raise ValueError(f"Archivo {name} no proporcionado.")
            
            content = file_obj.getbuffer()
            is_gz = file_obj.name.endswith(".gz")
            
            if is_gz:
                try:
                    content = gzip.decompress(content)
                except gzip.BadGzipFile:
                    raise ValueError(f"Archivo {name} ({file_obj.name}) parece .gz pero no se pudo descomprimir.")
            
            if expected_start: # Validación específica para MatrixMarket
                with io.BytesIO(content) as f_content:
                    first_line = f_content.readline().decode(errors='ignore')
                    if not first_line.startswith(expected_start):
                        raise ValueError(f"Archivo {name} ({file_obj.name}) no es válido. Esperaba '{expected_start}', obtuvo '{first_line[:50]}...'")
            else: # Validación genérica para TSV
                with io.BytesIO(content) as f_content:
                    # Leer algunas líneas para ver si es un TSV válido
                    # comment='%' puede ayudar si hay líneas de encabezado de MatrixMarket en features.tsv
                    pd.read_csv(f_content, sep="\t", header=None, nrows=5, comment='%', on_bad_lines='warn') 
        return True
    except ValueError as ve: # Errores de validación específicos
        st.error(f"Error de validación: {ve}")
        return False
    except Exception as e: # Otros errores de lectura/procesamiento
        st.error(f"Error procesando archivo para validación: {e}")
        # st.error(traceback.format_exc()) # Descomentar para depuración más profunda
        return False


def load_10x_data(matrix_file, features_file, barcodes_file, sample_name):
    """Carga datos 10x en un objeto AnnData."""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Guardar archivos con los nombres que Scanpy espera, manteniendo .gz si existe
        file_map = {
            "matrix.mtx": matrix_file,
            "features.tsv": features_file, # sc.read_10x_mtx también busca genes.tsv
            "barcodes.tsv": barcodes_file
        }
        for base_name, uploaded_file in file_map.items():
            filename_in_temp = base_name
            if uploaded_file.name.endswith(".gz"):
                filename_in_temp += ".gz" # Mantener la extensión .gz
            
            with open(os.path.join(temp_dir, filename_in_temp), "wb") as f:
                f.write(uploaded_file.getbuffer()) # Escribir el buffer original (comprimido o no)
        
        adata = sc.read_10x_mtx(temp_dir, var_names='gene_symbols', cache=False)
        adata.var_names_make_unique()
        adata.obs['sample'] = sample_name
        adata.X = sparse.csr_matrix(adata.X)
    return adata

def fig_to_bytes(fig, format='png'):
    """Convierte figura Matplotlib a bytes para descarga."""
    buf = io.BytesIO()
    fig.savefig(buf, format=format, bbox_inches='tight', dpi=300)
    buf.seek(0)
    return buf.getvalue()

def suggest_genes(input_genes, valid_genes_lower_map):
    """Sugiere genes similares (insensible a mayúsculas) si no se encuentran en el dataset."""
    suggestions = {}
    # valid_genes_original_case = list(valid_genes_lower_map.values()) # Nombres originales para sugerencias
    valid_genes_lower_list = list(valid_genes_lower_map.keys())     # Nombres en minúscula para búsqueda

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
    "leiden_flavor": "igraph", "umap_init_pos": "pca", "umap_n_neighbors": 15, "umap_min_dist": 0.5 # UMAP init pca
}
for key, value in default_values.items():
    if key not in st.session_state:
        st.session_state[key] = value

# --- Sidebar ---
with st.sidebar:
    st.image("https://www.biogenouest.org/wp-content/uploads/2020/02/Logo-ScanPy-195x150.png", width=100) # Ejemplo de logo
    st.header("Configuración del Análisis")
    
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
        
        leiden_flavors = ["igraph", "leidenalg"]
        st.session_state.leiden_flavor = st.selectbox("Backend Leiden", leiden_flavors, 
                                                      index=leiden_flavors.index(st.session_state.leiden_flavor) if st.session_state.leiden_flavor in leiden_flavors else 0, 
                                                      key="leiden_flavor_select", help="'igraph' es generalmente recomendado.")
        
        st.subheader("Parámetros UMAP (Avanzado)")
        umap_init_options = ["spectral", "random", "pca"]
        st.session_state.umap_init_pos = st.selectbox("Inicialización UMAP", umap_init_options, 
                                                      index=umap_init_options.index(st.session_state.umap_init_pos) if st.session_state.umap_init_pos in umap_init_options else 1, # Default a 'random'
                                                      key="umap_init_select", help="'random' puede ser más estable con algunas versiones de NumPy.")
        st.session_state.umap_n_neighbors = st.slider("Nº Vecinos UMAP (para embedding)", 2, 200, st.session_state.umap_n_neighbors, key="umap_n_neighbors_slider", help="Controla el balance entre estructura local y global en UMAP.")
        st.session_state.umap_min_dist = st.slider("Distancia Mínima UMAP", 0.0, 1.0, st.session_state.umap_min_dist, 0.01, key="umap_min_dist_slider", help="Controla cuán agrupados estarán los puntos en UMAP.")


    if st.session_state.analysis_done and st.session_state.adata_processed is not None:
        with st.expander("3. Análisis de Expresión Diferencial (DEA)"):
            # ... (Sección DEA, sin cambios mayores, pero asegurando que use adata_for_dea_config consistentemente) ...
            # ... y que la columna 'condition_temp_dea' se maneje bien.
            # ... (El código de DEA de tu versión anterior era bastante bueno, lo mantendré similar)
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
                # ... (Lógica de selectbox para dea_group1 y dea_group2, dea_cluster_scope, etc. como en tu última versión) ...
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

                # Crear columna temporal para DEA en una copia para evitar modificar el original innecesariamente
                adata_for_dea_preview = adata_for_dea_config.copy()
                adata_for_dea_preview.obs['condition_temp_dea'] = adata_for_dea_preview.obs['sample'].map(st.session_state.condition_assignments)
                
                # Filtrar NAs en condition_temp_dea antes del groupby para la vista previa
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
                                # Usar la copia con la columna 'condition_temp_dea' ya creada
                                adata_filtered_for_dea = adata_for_dea_preview[
                                    adata_for_dea_preview.obs['condition_temp_dea'].isin([st.session_state.dea_group1, st.session_state.dea_group2]) &
                                    adata_for_dea_preview.obs['condition_temp_dea'].notna() # Asegurar que no haya NAs
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


    # Botones principales de carga y ejecución del pipeline
    all_files_provided = True
    if st.session_state.num_samples > 0:
        for i in range(st.session_state.num_samples):
            if not (st.session_state.sample_files.get(f"matrix_file_{i}") and \
                    st.session_state.sample_files.get(f"features_file_{i}") and \
                    st.session_state.sample_files.get(f"barcodes_file_{i}")):
                all_files_provided = False
                break
    else:
        all_files_provided = False # No hay muestras para cargar

    if all_files_provided:
        if st.button("Cargar y Concatenar Datos", key="load_concat_btn", type="primary"):
            validation_messages = []
            all_valid_globally = True
            for i in range(st.session_state.num_samples):
                sample_name_val = st.session_state.sample_files[f"sample_name_{i}"]
                validation_messages.append(f"Validando Muestra {i+1} ({sample_name_val}):")
                if not validate_10x_files(
                    st.session_state.sample_files[f"matrix_file_{i}"],
                    st.session_state.sample_files[f"features_file_{i}"],
                    st.session_state.sample_files[f"barcodes_file_{i}"]
                ):
                    all_valid_globally = False
                    validation_messages.append(f" -> Archivos de Muestra {i+1} ({sample_name_val}) NO son válidos.")
                else:
                    validation_messages.append(f" -> Archivos de Muestra {i+1} ({sample_name_val}) VALIDADOS.")
            
            # Usar st.expander para mostrar mensajes de validación solo si hay errores o para confirmación
            with st.expander("Registro de Validación de Archivos", expanded=not all_valid_globally):
                for msg in validation_messages:
                    if "NO son válidos" in msg: st.error(msg)
                    else: st.info(msg)
            # Si todos los archivos son válidos, proceder a cargar y concatenar
            if all_valid_globally:
                with st.spinner("Cargando y concatenando datos..."):
                    try:
                        adatas_dict = {} # Usar un diccionario
                        sample_names_ordered = [] # Para mantener el orden si es necesario para otras cosas

                        for i in range(st.session_state.num_samples):
                            sample_name_user = st.session_state.sample_files[f"sample_name_{i}"]
                            sample_names_ordered.append(sample_name_user)

                            adata_sample = load_10x_data(
                                st.session_state.sample_files[f"matrix_file_{i}"],
                                st.session_state.sample_files[f"features_file_{i}"],
                                st.session_state.sample_files[f"barcodes_file_{i}"],
                                sample_name_user # load_10x_data ya añade .obs['sample']
                            )
                            if adata_sample is None:
                                raise ValueError(f"Fallo al cargar Muestra {i+1} ({sample_name_user}).")
                            
                            # No necesitas añadir .obs['sample'] aquí si load_10x_data lo hace,
                            # pero es crucial que la columna que concat usará para el label
                            # contenga estos nombres. O, mejor, usar el diccionario.
                            adatas_dict[sample_name_user] = adata_sample # Clave es el nombre de la muestra

                        # Concatenar usando el diccionario.
                        # El parámetro 'label' ahora nombrará la columna que contiene las CLAVES del diccionario.
                        st.session_state.adata_raw = ad.concat(
                            adatas_dict, # Pasar el diccionario
                            label='sample',    # La nueva columna se llamará 'sample' y contendrá los nombres de muestra
                            index_unique='-', 
                            join='outer', 
                            fill_value=0
                        )
                        # Ya no necesitas la línea:
                        # st.session_state.adata_raw.obs['sample'] = st.session_state.adata_raw.obs['sample_batch_id'].astype(str)
                        # porque .obs['sample'] ahora tiene los nombres correctos.

                        # Verificación (opcional, para DEBUG)
                        print("DEBUG: Valores únicos en st.session_state.adata_raw.obs['sample'] después de concat:", 
                              st.session_state.adata_raw.obs['sample'].unique())

                        st.session_state.adata_processed = None 
                        st.session_state.analysis_done = False
                        # ... (resetear otros estados) ...
                        st.success(f"Carga completada: {st.session_state.adata_raw.n_obs} células, {st.session_state.adata_raw.n_vars} genes.")
                    except Exception as e_load:
                        st.error(f"Error durante la carga: {e_load}")
                        st.error(traceback.format_exc())
                        st.session_state.adata_raw = None

                # Si la validación falla, no se carga nada
                # else: # El mensaje de error ya se mostró en el expander
            #     st.error("Algunos archivos no pasaron la validación. Por favor, revisa los mensajes de error detallados arriba.")

    elif st.session_state.num_samples > 0 :
        st.warning("Por favor, sube todos los archivos para cada muestra para habilitar la carga.")

    if st.session_state.adata_raw is not None:
        if st.button("Ejecutar Pipeline Principal", key="run_pipeline_btn", type="primary"):
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
                    progress_bar.progress(12) # Ajustado a 8 pasos

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
                    # n_neighbors debe ser < n_obs_hvg. Si n_obs_hvg es 2, max n_neighbors es 1.
                    valid_max_n_neighbors = n_obs_hvg - 1 if n_obs_hvg >1 else 1 
                    if current_n_neighbors_val > valid_max_n_neighbors: current_n_neighbors_val = valid_max_n_neighbors; st.warning(f"Nº Vecinos ajustado a {current_n_neighbors_val}")
                    if current_n_neighbors_val <= 0: current_n_neighbors_val = 1; st.warning(f"Nº Vecinos <=0, ajustado a {current_n_neighbors_val}")

                    print(f"DEBUG Neighbors: Usando n_neighbors={current_n_neighbors_val}, n_pcs={current_n_pcs_val}")
                    sc.pp.neighbors(st.session_state.adata_hvg_subset, n_neighbors=current_n_neighbors_val, n_pcs=current_n_pcs_val, random_state=0)
                    progress_bar.progress(75) # PCA y Neighbors hechos

                    # 7. UMAP
                    status_text.text("Paso 7/8: Cálculo de UMAP...")
                    print(f"DEBUG UMAP: init={st.session_state.umap_init_pos}, n_neighbors={st.session_state.umap_n_neighbors}, min_dist={st.session_state.umap_min_dist}")
                    try:
                        import umap # Mover import aquí por si no se usa siempre
                        reducer = umap.UMAP(
                            n_neighbors=int(st.session_state.umap_n_neighbors), n_components=2, metric='euclidean',
                            min_dist=float(st.session_state.umap_min_dist), init=st.session_state.umap_init_pos,
                            random_state=42, verbose=False 
                        )
                        X_pca_for_umap = st.session_state.adata_hvg_subset.obsm['X_pca']
                        embedding = reducer.fit_transform(X_pca_for_umap)
                        st.session_state.adata_hvg_subset.obsm['X_umap'] = embedding
                        print("DEBUG UMAP: API directa completada.")
                    except Exception as e_umap_direct_api:
                        st.warning(f"UMAP con API directa falló ({e_umap_direct_api}). Intentando con sc.tl.umap como fallback...")
                        try:
                            umap_init_fallback = st.session_state.umap_init_pos if st.session_state.umap_init_pos != 'pca' else 'spectral'
                            sc.tl.umap(st.session_state.adata_hvg_subset, 
                                       n_neighbors=int(st.session_state.umap_n_neighbors), 
                                       min_dist=float(st.session_state.umap_min_dist),
                                       init_pos=umap_init_fallback, 
                                       random_state=0) 
                            print("DEBUG UMAP: sc.tl.umap (fallback) completado.")
                        except Exception as e_sc_umap_fallback:
                             st.error(f"Fallback sc.tl.umap también falló: {e_sc_umap_fallback}. UMAP no se calculará.")
                             st.session_state.adata_hvg_subset.obsm.pop('X_umap', None)
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
                    
                    # Marcadores
                    adata.obs['leiden_clusters'] = adata.obs['leiden_clusters'].astype('category')
                    sc.tl.rank_genes_groups(adata, 'leiden_clusters', method='wilcoxon', key_added='rank_genes_leiden_clusters', use_raw=False)
                    res_markers = adata.uns.get('rank_genes_leiden_clusters', {}) 
                    marker_data_list = []
                    
                    # Verificar que res_markers no esté vacío y tenga los campos esperados
                    # El .dtype.names ya nos da los grupos que tienen 'names'
                    if res_markers and 'names' in res_markers and hasattr(res_markers['names'], 'dtype') and res_markers['names'].dtype.names is not None:
                        cluster_ids_with_names = res_markers['names'].dtype.names
                        
                        for grp_marker in cluster_ids_with_names:
                            # Verificar que este grp_marker (ID de clúster) también sea un campo válido
                            # en los otros arrays estructurados.
                            if ('scores' in res_markers and hasattr(res_markers['scores'], 'dtype') and grp_marker in res_markers['scores'].dtype.names and \
                                'logfoldchanges' in res_markers and hasattr(res_markers['logfoldchanges'], 'dtype') and grp_marker in res_markers['logfoldchanges'].dtype.names and \
                                'pvals_adj' in res_markers and hasattr(res_markers['pvals_adj'], 'dtype') and grp_marker in res_markers['pvals_adj'].dtype.names):
                                
                                num_avail_markers = len(res_markers['names'][grp_marker])
                                markers_to_fetch = min(st.session_state.n_top_markers, num_avail_markers)
                                
                                for i_marker in range(markers_to_fetch):
                                    # Comprobación adicional por si alguno de los arrays es más corto inesperadamente para este índice
                                    if i_marker < len(res_markers['scores'][grp_marker]) and \
                                       i_marker < len(res_markers['logfoldchanges'][grp_marker]) and \
                                       i_marker < len(res_markers['pvals_adj'][grp_marker]):
                                        marker_data_list.append({
                                            'Cluster': grp_marker, 
                                            'Rank': i_marker + 1, 
                                            'Gene': res_markers['names'][grp_marker][i_marker], 
                                            'Score': res_markers['scores'][grp_marker][i_marker],
                                            'Log2FC': res_markers['logfoldchanges'][grp_marker][i_marker], 
                                            'P-Value Adj': res_markers['pvals_adj'][grp_marker][i_marker]
                                        })
                                    else:
                                        print(f"Warning: Discrepancia en la longitud de arrays de marcadores para el clúster '{grp_marker}' en el índice {i_marker}.")
                            else:
                                print(f"Warning: Faltan campos de datos ('scores', 'logfoldchanges' o 'pvals_adj') para el grupo de marcadores '{grp_marker}' o el grupo no es un campo en todos ellos.")
                    else:
                        print("Warning: 'rank_genes_leiden_clusters' no tiene la estructura esperada en .uns o no tiene 'names'.")

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
                        # Resetear estados de resultados
                        for key_to_reset in ["marker_genes_df", "dea_results_df", "dea_comparison_str", 
                                             "dea_group1", "dea_group2", "dea_cluster_scope",
                                             "dea_n_genes_display", "dea_lfc_cutoff", "dea_pval_cutoff"]: # Incluir todos los relevantes
                            if key_to_reset in default_values: 
                                st.session_state[key_to_reset] = default_values[key_to_reset]
                            elif key_to_reset in st.session_state: # Si no está en defaults pero existe, poner a None
                                 st.session_state[key_to_reset] = None


# --- Sección de Resultados ---
# (El código de esta sección se mantiene igual que en tu última versión, con las mejoras para
# el explorador de genes y el manejo de 'X_umap' potencialmente ausente.
# Lo omito aquí por brevedad, pero debes pegarlo completo.)
st.header("Resultados del Análisis")
st.subheader("🔬 Explorador de Expresión Génica")
st.session_state.gene_explorer_input = st.text_area(
    "Ingresa nombres de genes (separados por coma, espacio o nueva línea):", 
    value=st.session_state.gene_explorer_input, 
    key="gene_explorer_main_input",
    height=100
)

if st.session_state.analysis_done and st.session_state.adata_processed is not None:
    adata_display = st.session_state.adata_processed 
    
    valid_genes_lower_map_display = {g.lower(): g for g in adata_display.var_names}

    try:
        # Escribir en un archivo temporal
        with tempfile.NamedTemporaryFile(delete=False, suffix=".h5ad") as tmp_h5ad_file:
            filepath_h5ad = tmp_h5ad_file.name
        # El archivo se cierra automáticamente al salir del with
        adata_display.write_h5ad(filepath_h5ad) # Escribe en el path obtenido

        with open(filepath_h5ad, "rb") as f_h5ad_read:
            st.sidebar.download_button(
                "Descargar AnnData Procesado (.h5ad)", 
                f_h5ad_read.read(), 
                f"processed_adata_scRNAseq_{pd.Timestamp.now().strftime('%Y%m%d_%H%M')}.h5ad",
                "application/octet-stream",
                key="download_adata_button_v3" # Nueva key para evitar conflictos
            )
        os.remove(filepath_h5ad) 
    except Exception as e_dl:
        st.sidebar.error(f"Error al preparar descarga de AnnData: {e_dl}")
        # Considerar si el archivo temporal debe eliminarse aquí también si existe
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

    tab_titles_main = ["📊 UMAPs", "🔬 Marcadores de Clúster", "🧬 Plots de QC", "📈 Análisis Diferencial", "🧬 Explorador de Genes", "ℹ️ Info del Dataset"]
    tab_umaps_main, tab_markers_main, tab_qc_main, tab_dea_main, tab_gene_explorer_main, tab_info_main = st.tabs(tab_titles_main)
    
    with tab_umaps_main:
        if 'X_umap' not in adata_display.obsm:
            st.warning("Los datos de UMAP no están disponibles. El cálculo de UMAP pudo haber fallado o no se completó.")
        else:
            st.subheader("UMAP coloreado por Clústeres de Leiden")
            if 'leiden_clusters' in adata_display.obs:
                fig_umap_clusters, ax_clusters = plt.subplots(figsize=(7, 6))
                sc.pl.umap(adata_display, color='leiden_clusters', legend_loc='on data', ax=ax_clusters, show=False, title=f"Clusters Leiden (Res: {st.session_state.leiden_res})")
                st.pyplot(fig_umap_clusters)
                st.download_button("Descargar UMAP Clústeres (PNG)", fig_to_bytes(fig_umap_clusters), "umap_leiden_clusters.png", "image/png", key="dl_umap_c_v2")
                plt.close(fig_umap_clusters)
            else: st.warning("Clusters de Leiden no disponibles para el plot UMAP.")

            if 'sample' in adata_display.obs:
                st.subheader("UMAP coloreado por Muestra")
                fig_umap_sample, ax_sample = plt.subplots(figsize=(7, 6))
                sc.pl.umap(adata_display, color='sample', ax=ax_sample, show=False, title="Por Muestra")
                st.pyplot(fig_umap_sample)
                st.download_button("Descargar UMAP Muestra (PNG)", fig_to_bytes(fig_umap_sample), "umap_by_sample.png", "image/png", key="dl_umap_s_v2")
                plt.close(fig_umap_sample)

            st.subheader("UMAPs por Muestra (Facetado, Coloreado por Clúster)")
            if 'sample' in adata_display.obs and 'leiden_clusters' in adata_display.obs:
                try:
                    unique_samples_facet = sorted(adata_display.obs['sample'].astype('category').cat.categories.tolist())
                    n_samples_facet = len(unique_samples_facet)
                    if n_samples_facet > 0:
                        cols_facet = min(n_samples_facet, 3); rows_facet = (n_samples_facet + cols_facet - 1) // cols_facet
                        fig_facet, axes_facet = plt.subplots(rows_facet, cols_facet, figsize=(cols_facet * 5, rows_facet * 4.5), squeeze=False)
                        axes_flat_facet = axes_facet.flatten()
                        idx_facet = 0 # Para asegurar que idx_facet esté definido antes del bucle de limpieza
                        for idx_facet, sample_val_facet in enumerate(unique_samples_facet):
                            if idx_facet < len(axes_flat_facet):
                                ax_curr_facet = axes_flat_facet[idx_facet]
                                adata_subset_facet = adata_display[adata_display.obs['sample'] == sample_val_facet].copy() # Copia para evitar warnings
                                if not adata_subset_facet.obs.empty:
                                    sc.pl.umap(adata_subset_facet, color='leiden_clusters', ax=ax_curr_facet, show=False, 
                                               title=f"Muestra: {sample_val_facet}", legend_loc='on data' if idx_facet == 0 else None, legend_fontsize=6)
                                else:
                                    ax_curr_facet.text(0.5, 0.5, f"M: {sample_val_facet}\n(Sin células)", ha='center', va='center'); ax_curr_facet.set_xticks([]); ax_curr_facet.set_yticks([])
                        for j_ax_empty in range(idx_facet + 1, len(axes_flat_facet)): fig_facet.delaxes(axes_flat_facet[j_ax_empty])
                        plt.tight_layout(); st.pyplot(fig_facet)
                        st.download_button("Descargar UMAPs Facetados (PNG)", fig_to_bytes(fig_facet), "umaps_faceted_by_sample.png", "image/png", key="dl_umaps_facet_v2")
                        plt.close(fig_facet)
                except Exception as e_facet: st.error(f"Error UMAPs facetados: {e_facet}")

    with tab_markers_main:
        if st.session_state.marker_genes_df is not None and not st.session_state.marker_genes_df.empty:
            st.subheader(f"Top {st.session_state.n_top_markers} Genes Marcadores por Clúster")
            st.dataframe(st.session_state.marker_genes_df)
            st.download_button(
                "Descargar Marcadores (CSV)", st.session_state.marker_genes_df.to_csv(index=False).encode('utf-8'),
                "cluster_markers.csv", "text/csv", key="dl_markers_v2"
            )
            
            st.subheader("Dot Plot de Genes Marcadores")
            df_markers_for_plot = st.session_state.marker_genes_df
            top_n_per_cluster_dotplot = st.slider("Nº top marcadores por clúster para Dot Plot", 1, 5, 2, key="dotplot_n_markers_slider") # Permitir al usuario elegir
            genes_for_marker_dotplot_list = []
            if 'Cluster' in df_markers_for_plot.columns and 'Rank' in df_markers_for_plot.columns and 'Gene' in df_markers_for_plot.columns:
                for cluster_id_dotplot in sorted(df_markers_for_plot['Cluster'].astype(str).unique()): # Asegurar que Cluster ID es string
                    cluster_specific_markers = df_markers_for_plot[df_markers_for_plot['Cluster'] == cluster_id_dotplot]
                    genes_for_marker_dotplot_list.extend(cluster_specific_markers.nsmallest(top_n_per_cluster_dotplot, 'Rank')['Gene'].tolist())
                unique_genes_for_dotplot = list(dict.fromkeys(genes_for_marker_dotplot_list))

                if unique_genes_for_dotplot and 'leiden_clusters' in adata_display.obs:
                    try:
                        num_clusters_dotplot = adata_display.obs['leiden_clusters'].nunique()
                        fig_dot_markers, ax_dot_markers = plt.subplots(figsize=(max(8, len(unique_genes_for_dotplot) * 0.6), max(5, num_clusters_dotplot * 0.45)))
                        sc.pl.dotplot(adata_display, unique_genes_for_dotplot, groupby='leiden_clusters', ax=ax_dot_markers, show=False, standard_scale='var', use_raw=False)
                        plt.xticks(rotation=90); plt.tight_layout(); st.pyplot(fig_dot_markers)
                        st.download_button("Descargar Dot Plot Marcadores (PNG)", fig_to_bytes(fig_dot_markers), "dotplot_cluster_markers.png", "image/png", key="dl_dotplot_markers_v2")
                        plt.close(fig_dot_markers)
                    except Exception as e_dot_m: st.error(f"Error dot plot marcadores: {e_dot_m}")
                else: st.info("No hay genes marcadores para el dot plot o faltan clusters.")
            else: st.info("Columnas ('Cluster', 'Rank', 'Gene') no encontradas para dot plot de marcadores.")
        else: st.info("No se han calculado genes marcadores.")


    with tab_qc_main:
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
                        st.download_button(f"Descargar {metric_title.replace(' ', '_')} (PNG)", fig_to_bytes(fig_qc_violin), f"qc_violin_{metric_key}.png", "image/png", key=f"dl_qc_vln_{metric_key}_v2")
                        plt.close(fig_qc_violin)
                    except Exception as e_qc_plot: st.error(f"Error violin QC para {metric_title}: {e_qc_plot}")
                else: st.warning(f"Métrica QC '{metric_key}' no encontrada.")
        else: st.warning("Columna 'sample' no encontrada para plots de QC.")

    with tab_dea_main:
        # ... (Código de la pestaña DEA como en tu última versión, era bastante bueno) ...
        # ... asegurando que usa adata_display para cualquier plot o dataframe basado en el AnnData procesado.
        # ... y que los nombres de archivo para descarga sean robustos.
        st.subheader("Resultados del Análisis de Expresión Diferencial")
        if st.session_state.dea_results_df is not None and not st.session_state.dea_results_df.empty:
            st.markdown(f"**Comparación Actual:** `{st.session_state.dea_comparison_str}`")
            st.dataframe(st.session_state.dea_results_df.head(st.session_state.dea_n_genes_display))
            
            csv_dea_filename_safe = "".join(c if c.isalnum() or c in (' ', '_', '-') else '_' for c in st.session_state.dea_comparison_str).rstrip()
            csv_dea_filename = f"dea_results_{csv_dea_filename_safe.replace(' vs ','_vs_').replace(' ','_')}.csv"
            st.download_button(
                "Descargar Tabla DEA Completa (CSV)", 
                st.session_state.dea_results_df.to_csv(index=False).encode('utf-8'),
                csv_dea_filename, "text/csv", key="dl_dea_table_csv_v2"
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
                fig_px_volcano.add_hline(y=-np.log10(st.session_state.dea_pval_cutoff), line_dash="dash", line_color="black", annotation_text=f"P.adj={st.session_state.dea_pval_cutoff}")
                fig_px_volcano.add_vline(x=st.session_state.dea_lfc_cutoff, line_dash="dash", line_color="black", annotation_text=f"LFC={st.session_state.dea_lfc_cutoff}")
                fig_px_volcano.add_vline(x=-st.session_state.dea_lfc_cutoff, line_dash="dash", line_color="black", annotation_text=f"LFC={-st.session_state.dea_lfc_cutoff}")
                st.plotly_chart(fig_px_volcano, use_container_width=True)
                
                html_volcano_buffer = io.StringIO()
                fig_px_volcano.write_html(html_volcano_buffer)
                html_volcano_filename_safe = "".join(c if c.isalnum() or c in (' ', '_', '-') else '_' for c in st.session_state.dea_comparison_str).rstrip()
                html_volcano_filename = f"volcano_plot_{html_volcano_filename_safe.replace(' vs ','_vs_').replace(' ','_')}.html"
                st.download_button("Descargar Volcano Plot (HTML)", html_volcano_buffer.getvalue(), html_volcano_filename, "text/html", key="dl_volcano_html_v2")
            except Exception as e_volc:
                st.error(f"Error generando Volcano Plot: {e_volc}")
                st.error(traceback.format_exc())
        elif st.session_state.analysis_done:
            st.info("No hay resultados de DEA para mostrar. Ejecuta el Análisis Diferencial.")


    with tab_gene_explorer_main:
        # ... (código del explorador de genes como estaba, asegurando que usa adata_display) ...
        # ... y que maneja si X_umap no existe.
        st.subheader("Visualización de Expresión para Genes Específicos")
        if not genes_to_visualize_list:
            st.info("Ingresa nombres de genes válidos en el campo de texto de arriba para visualizarlos.")
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
                    idx_ge_plot = 0 # Para asegurar que idx_ge_plot esté definido
                    for idx_ge_plot, gene_name_plot in enumerate(genes_to_visualize_list):
                        if idx_ge_plot < len(axes_flat_ge_umaps):
                            ax_ge_curr = axes_flat_ge_umaps[idx_ge_plot]
                            try:
                                sc.pl.umap(adata_display, color=gene_name_plot, ax=ax_ge_curr, show=False, title=gene_name_plot, cmap='viridis', use_raw=False)
                            except Exception as e_ge_umap_plot: 
                                ax_ge_curr.text(0.5, 0.5, f"Error plot\n{gene_name_plot}", ha='center', va='center', color='red'); ax_ge_curr.set_xticks([]); ax_ge_curr.set_yticks([])
                                print(f"Error ploteando UMAP para gen {gene_name_plot}: {e_ge_umap_plot}") # Log a consola
                    for j_ge_empty_ax in range(idx_ge_plot + 1, len(axes_flat_ge_umaps)): fig_ge_umaps.delaxes(axes_flat_ge_umaps[j_ge_empty_ax])
                    plt.tight_layout(); st.pyplot(fig_ge_umaps)
                    st.download_button("Descargar UMAPs de Genes (PNG)", fig_to_bytes(fig_ge_umaps), "gene_explorer_umaps.png", "image/png", key="dl_ge_umaps_png_v2")
                    plt.close(fig_ge_umaps)
            
            # (Resto de plots del explorador de genes: violines, dotplot...)
            if 'leiden_clusters' in adata_display.obs and genes_to_visualize_list:
                st.markdown("#### Diagramas de Violín por Clúster de Leiden")
                try:
                    genes_for_violin_cluster = genes_to_visualize_list[:min(5, len(genes_to_visualize_list))]
                    n_clusters_violin = adata_display.obs['leiden_clusters'].nunique()
                    fig_ge_violins_cl, ax_ge_violins_cl = plt.subplots(figsize=(max(7, n_clusters_violin * 0.8), 5))
                    sc.pl.violin(adata_display, keys=genes_for_violin_cluster, groupby='leiden_clusters', rotation=45, ax=ax_ge_violins_cl, show=False, use_raw=False, cut=0, multi_panel=len(genes_for_violin_cluster)>1)
                    plt.tight_layout(); st.pyplot(fig_ge_violins_cl)
                    st.download_button("Violines por Clúster (PNG)", fig_to_bytes(fig_ge_violins_cl), "ge_violins_cluster.png", key="dl_ge_vln_cl_v2")
                    plt.close(fig_ge_violins_cl)
                except Exception as e_ge_vln_cl: st.error(f"Error violines por clúster: {e_ge_vln_cl}")
            
            if 'condition_temp_dea' in adata_display.obs and adata_display.obs['condition_temp_dea'].nunique() > 1 and genes_to_visualize_list:
                st.markdown("#### Diagramas de Violín por Condición (definida en DEA)")
                # ... (código para violines por condición)
                pass

            if len(genes_to_visualize_list) > 0 and 'leiden_clusters' in adata_display.obs:
                st.markdown("#### Dot Plot de Genes Seleccionados por Clúster")
                # ... (código para dotplot de genes del explorador)
                pass


    with tab_info_main:
        # ... (código de info como estaba, usando adata_display) ...
        st.subheader("Información del Dataset Procesado")
        st.write(f"Total de Células (post-QC): {adata_display.n_obs}")
        st.write(f"Total de Genes (post-QC): {adata_display.n_vars}")
        if st.session_state.adata_hvg_subset is not None: # Corregido para reflejar el nombre correcto
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
    # Mensajes de bienvenida o estado inicial
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
st.sidebar.info("Analizador scRNA-seq v0.8. Basado en Scanpy.")
st.sidebar.markdown("Se recomienda crear un entorno virtual con versiones compatibles de las bibliotecas (ej: `numpy<2.0` si se experimentan errores con UMAP).")
st.sidebar.markdown("Si tienes problemas, consulta la [documentación de Scanpy](https://scanpy.readthedocs.io/en/stable/) o el [repositorio de GitHub]).")
st.sidebar.markdown("Para más información, visita el [repositorio de GitHub]).")
st.sidebar.markdown("**Desarrollado por:** Pedro Botías - pbotias@ucm.es - https://github.com/pbotiast/scRNASeq")
st.sidebar.markdown("**Licencia:** Licencia MIT - https://opensource.org/licenses/MIT")
st.sidebar.markdown("**Fecha:** [05/05/2025] - [12/05/2025]")  
st.sidebar.markdown("**Versión:** 0.8")
st.sidebar.markdown("**Última Actualización:** 2025-05-12")
st.sidebar.markdown("**Notas:** Esta aplicación es un prototipo y puede contener errores. Usa bajo tu propio riesgo.")
st.sidebar.markdown("**Disclaimer:** Esta aplicación es un prototipo y puede contener errores. Usa bajo tu propio riesgo.")



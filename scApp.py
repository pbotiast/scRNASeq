# -*- coding: utf-8 -*-
"""
Streamlit app para análisis de datos de Single-Cell RNA-seq con múltiples muestras.
Este script permite cargar datos de múltiples muestras, realizar análisis de calidad, normalización, PCA, UMAP,
clustering y análisis diferencial.
"""
import streamlit as st
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import pandas as pd
import os
import tempfile
import io
import traceback
import plotly.express as px

# Configuración de la página de Streamlit
st.set_page_config(layout="wide")
st.title("Analizador Interactivo de Single-Cell RNA-seq")

# --- Funciones Auxiliares ---
def load_10x_data(matrix_file, features_file, barcodes_file, sample_name_for_adata):
    with tempfile.TemporaryDirectory() as temp_dir:
        with open(os.path.join(temp_dir, matrix_file.name), "wb") as f: f.write(matrix_file.getbuffer())
        with open(os.path.join(temp_dir, features_file.name), "wb") as f: f.write(features_file.getbuffer())
        with open(os.path.join(temp_dir, barcodes_file.name), "wb") as f: f.write(barcodes_file.getbuffer())
        adata_sample = sc.read_10x_mtx(temp_dir, var_names='gene_symbols', cache=True)
        adata_sample.var_names_make_unique()
        adata_sample.obs['sample'] = sample_name_for_adata
    return adata_sample

def fig_to_bytes(fig, format='png'):
    buf = io.BytesIO()
    fig.savefig(buf, format=format, bbox_inches='tight', dpi=300); buf.seek(0)
    return buf.getvalue()

# --- Inicialización Centralizada de st.session_state ---
default_values = {
    "num_samples_value": 1, "sample_files": {},
    "min_genes_val": 200, "min_cells_val": 3, "mito_prefix_val": "MT-", "max_mito_pct_val": 10,
    "n_top_genes_hvg_val": 2000, "n_pcs_val": 30, "leiden_res_val": 0.8, "n_top_markers_val": 5,
    "adata_combined_raw": None, "adata_processed": None, "adata_hvg_filtered_intermediate": None,
    "analysis_done": False, "marker_genes_df": None,
    "condition_assignments": {}, "dea_group1": None, "dea_group2": None,
    "dea_cluster_scope": "Todos los Clústeres", "dea_n_genes_display": 25,
    "dea_lfc_cutoff": 0.5, "dea_pval_cutoff": 0.05,
    "dea_results_df": None, "dea_comparison_str": "",
    "gene_explorer_input": ""
}
for key, value in default_values.items():
    if key not in st.session_state:
        st.session_state[key] = value

# --- Sidebar ---
with st.sidebar.expander("1. Carga de Datos", expanded=True):
    st.session_state.num_samples_value = st.number_input(
        "Número de muestras a cargar", min_value=1, max_value=10,
        value=st.session_state.num_samples_value, step=1, key="num_samples_widget_main"
    )
    for i in range(st.session_state.num_samples_value):
        st.subheader(f"Muestra {i+1}")
        s_name_key, m_key, f_key, b_key = f"sample_name_{i}", f"matrix_file_{i}", f"features_file_{i}", f"barcodes_file_{i}"
        
        # Asegurar que las sub-claves existen en sample_files
        if s_name_key not in st.session_state.sample_files: st.session_state.sample_files[s_name_key] = f"Muestra{i+1}"
        if m_key not in st.session_state.sample_files: st.session_state.sample_files[m_key] = None
        if f_key not in st.session_state.sample_files: st.session_state.sample_files[f_key] = None
        if b_key not in st.session_state.sample_files: st.session_state.sample_files[b_key] = None

        st.session_state.sample_files[s_name_key] = st.text_input(f"Nombre Muestra {i+1}", value=st.session_state.sample_files[s_name_key], key=s_name_key + "_exp_in")
        st.session_state.sample_files[m_key] = st.file_uploader(f"Matrix.mtx (M{i+1})", type=["mtx", "gz"], key=m_key + "_exp_in")
        st.session_state.sample_files[f_key] = st.file_uploader(f"Features.tsv (M{i+1})", type=["tsv", "gz"], key=f_key + "_exp_in")
        st.session_state.sample_files[b_key] = st.file_uploader(f"Barcodes.tsv (M{i+1})", type=["tsv", "gz"], key=b_key + "_exp_in")

with st.sidebar.expander("2. Parámetros de Pipeline Principal"):
    st.session_state.min_genes_val = st.slider("Mínimo genes/célula", 50, 1000, st.session_state.min_genes_val, help="Filtra células con menos genes detectados.")
    st.session_state.min_cells_val = st.slider("Mínimo células/gen", 1, 50, st.session_state.min_cells_val, help="Filtra genes expresados en menos células.")
    st.session_state.mito_prefix_val = st.text_input("Prefijo genes mitocondriales", st.session_state.mito_prefix_val)
    st.session_state.max_mito_pct_val = st.slider("Máx % mito", 1, 100, st.session_state.max_mito_pct_val, help="Filtra células con alto porcentaje de cuentas mitocondriales.")
    st.session_state.n_top_genes_hvg_val = st.slider("Nº HVGs", 500, 5000, st.session_state.n_top_genes_hvg_val, help="Número de genes altamente variables a seleccionar.")
    st.session_state.n_pcs_val = st.slider("Nº PCs", 10, 100, st.session_state.n_pcs_val, help="Número de componentes principales para PCA y vecinos.")
    st.session_state.leiden_res_val = st.slider("Resolución Leiden", 0.1, 2.0, st.session_state.leiden_res_val, 0.1, help="Resolución para el clustering de Leiden.")
    st.session_state.n_top_markers_val = st.slider("Nº marcadores/clúster", 1, 20, st.session_state.n_top_markers_val, help="Número de genes marcadores a mostrar por clúster.")

if st.session_state.analysis_done and st.session_state.adata_processed is not None:
    with st.sidebar.expander("3. Análisis de Expresión Diferencial (DEA)"):
        adata_for_dea_setup = st.session_state.adata_processed
        available_samples_for_dea = sorted(adata_for_dea_setup.obs['sample'].unique().tolist())

        st.subheader("Asignar Muestras a Condiciones")
        current_assignments = st.session_state.condition_assignments.copy() # Trabajar con copia para actualizar
        for sample_name_dea in available_samples_for_dea:
            default_cond = current_assignments.get(sample_name_dea, f"Cond_{sample_name_dea}")
            current_assignments[sample_name_dea] = st.text_input(f"Condición para {sample_name_dea}", value=default_cond, key=f"cond_for_{sample_name_dea}_exp")
        st.session_state.condition_assignments = current_assignments # Actualizar el estado

        unique_conditions = sorted(list(set(st.session_state.condition_assignments.values())))
        valid_conditions_dea = [cond for cond in unique_conditions if cond and cond.strip()]

        if len(valid_conditions_dea) >= 2:
            st.subheader("Seleccionar Grupos para Comparación")
            col1_dea_exp, col2_dea_exp = st.columns(2)
            with col1_dea_exp:
                st.session_state.dea_group1 = st.selectbox("Grupo 1 (Referencia)", options=valid_conditions_dea, 
                                                           index=valid_conditions_dea.index(st.session_state.dea_group1) if st.session_state.dea_group1 in valid_conditions_dea else 0, 
                                                           key="dea_g1_select_exp_v2")
            with col2_dea_exp:
                options_g2_dea = [cond for cond in valid_conditions_dea if cond != st.session_state.dea_group1]
                if not options_g2_dea and valid_conditions_dea: options_g2_dea = valid_conditions_dea
                st.session_state.dea_group2 = st.selectbox("Grupo 2 (Comparación)", options=options_g2_dea, 
                                                           index=options_g2_dea.index(st.session_state.dea_group2) if st.session_state.dea_group2 in options_g2_dea and options_g2_dea else 0, 
                                                           key="dea_g2_select_exp_v2")

            available_clusters_dea = ["Todos los Clústeres"] + sorted(adata_for_dea_setup.obs['leiden_clusters'].unique().tolist())
            st.session_state.dea_cluster_scope = st.selectbox("Realizar DEA en:", options=available_clusters_dea, 
                                                              index=available_clusters_dea.index(st.session_state.dea_cluster_scope) if st.session_state.dea_cluster_scope in available_clusters_dea else 0,
                                                              key="dea_cluster_scope_select_exp_v2")
            
            st.session_state.dea_n_genes_display = st.slider("Nº genes DEA", 10, 200, st.session_state.dea_n_genes_display, key="dea_n_genes_slider_exp_v2")
            st.session_state.dea_lfc_cutoff = st.number_input("Log2FC cutoff", 0.0, value=st.session_state.dea_lfc_cutoff, step=0.1, key="dea_lfc_cutoff_exp_v2")
            st.session_state.dea_pval_cutoff = st.number_input("P-adj cutoff", 0.0, 1.0, value=st.session_state.dea_pval_cutoff, step=0.01, format="%.3f", key="dea_pval_cutoff_exp_v2")

            if st.button("Ejecutar Análisis Diferencial", key="run_dea_button_exp_v2"):
                if st.session_state.dea_group1 == st.session_state.dea_group2:
                    st.error("Grupos de comparación deben ser diferentes.")
                elif not st.session_state.dea_group1 or not st.session_state.dea_group2:
                    st.error("Seleccione ambos grupos para comparación.")
                else:
                    with st.spinner("Ejecutando Análisis Diferencial..."):
                        try:
                            adata_dea_run = st.session_state.adata_processed.copy()
                            adata_dea_run.obs['condition'] = adata_dea_run.obs['sample'].map(st.session_state.condition_assignments).astype('category')
                            cells_g1 = adata_dea_run.obs['condition'] == st.session_state.dea_group1
                            cells_g2 = adata_dea_run.obs['condition'] == st.session_state.dea_group2
                            adata_dea_subset_cond = adata_dea_run[cells_g1 | cells_g2].copy()
                            
                            target_adata_for_dea_run = adata_dea_subset_cond
                            scope_msg = ""
                            if st.session_state.dea_cluster_scope != "Todos los Clústeres":
                                s_scope = st.session_state.dea_cluster_scope; scope_msg = f" (Clúster {s_scope})"
                                if s_scope in target_adata_for_dea_run.obs['leiden_clusters'].unique():
                                    target_adata_for_dea_run = target_adata_for_dea_run[target_adata_for_dea_run.obs['leiden_clusters'] == s_scope].copy()
                                else: raise ValueError(f"Clúster {s_scope} no válido para condiciones dadas.")
                            
                            counts_by_cond = target_adata_for_dea_run.obs['condition'].value_counts()
                            if target_adata_for_dea_run.n_obs == 0 or len(counts_by_cond) < 2 or \
                               counts_by_cond.get(st.session_state.dea_group1,0) < 3 or counts_by_cond.get(st.session_state.dea_group2,0) < 3:
                                st.error(f"Insuficientes células/grupos para DEA{scope_msg}. Counts: {counts_by_cond.to_dict()}")
                                st.session_state.dea_results_df = None
                            else:
                                sc.tl.rank_genes_groups(target_adata_for_dea_run, 'condition', groups=[st.session_state.dea_group2], reference=st.session_state.dea_group1, method='wilcoxon', key_added='rank_genes_dea', use_raw=False)
                                res_dea = target_adata_for_dea_run.uns['rank_genes_dea']
                                g_name_uns = st.session_state.dea_group2
                                if g_name_uns in res_dea['names'].dtype.names:
                                    df_res_dea = pd.DataFrame({'Gene': res_dea['names'][g_name_uns], 'Log2FC': res_dea['logfoldchanges'][g_name_uns], 'P-Value': res_dea['pvals'][g_name_uns], 'P-Value Adj': res_dea['pvals_adj'][g_name_uns], 'Scores': res_dea['scores'][g_name_uns]})
                                    st.session_state.dea_results_df = df_res_dea.sort_values(by='P-Value Adj')
                                    st.session_state.dea_comparison_str = f"{st.session_state.dea_group2} vs {st.session_state.dea_group1}{scope_msg}"
                                    st.success(f"DEA{scope_msg} completado.")
                                else: st.error(f"Grupo '{g_name_uns}' no en resultados DEA."); st.session_state.dea_results_df = None
                        except Exception as e_run_dea: st.error(f"Error en DEA: {e_run_dea}"); st.error(traceback.format_exc()); st.session_state.dea_results_df = None
        elif st.session_state.analysis_done: st.info("Define al menos dos condiciones válidas para comparar.")

st.sidebar.markdown("---")
# Botones de acción principales
all_files_uploaded_check = True
if st.session_state.num_samples_value > 0:
    for i in range(st.session_state.num_samples_value):
        if not (st.session_state.sample_files.get(f"matrix_file_{i}") and st.session_state.sample_files.get(f"features_file_{i}") and st.session_state.sample_files.get(f"barcodes_file_{i}")):
            all_files_uploaded_check = False; break
else: all_files_uploaded_check = False

if all_files_uploaded_check:
    if st.sidebar.button("Cargar y Concatenar Datos", key="load_data_main_btn"):
        with st.spinner("Cargando y concatenando..."):
            try:
                adatas_list_main = [load_10x_data(st.session_state.sample_files[f"matrix_file_{i}"], st.session_state.sample_files[f"features_file_{i}"], st.session_state.sample_files[f"barcodes_file_{i}"], st.session_state.sample_files[f"sample_name_{i}"]) for i in range(st.session_state.num_samples_value)]
                st.session_state.adata_combined_raw = ad.concat(adatas_list_main, label='sample_batch_id', index_unique='-', join='outer', fill_value=0)
                st.session_state.adata_processed = None; st.session_state.analysis_done = False; st.session_state.marker_genes_df = None; st.session_state.dea_results_df = None
                st.success(f"Cargado: {st.session_state.adata_combined_raw.n_obs} células, {st.session_state.adata_combined_raw.n_vars} genes.")
            except Exception as e_load_main: st.error(f"Error carga: {e_load_main}"); st.error(traceback.format_exc()); st.session_state.adata_combined_raw = None
elif st.session_state.num_samples_value > 0: st.sidebar.warning("Sube todos los archivos de muestra.")

if st.session_state.adata_combined_raw is not None:
    if st.sidebar.button("Ejecutar Pipeline Principal", key="run_analysis_main_btn"):
        adata_pipeline = st.session_state.adata_combined_raw.copy()
        st.session_state.adata_hvg_filtered_intermediate = None
        with st.spinner("Ejecutando pipeline..."):
            progress_bar_main = st.progress(0); status_text_main = st.empty()
            try:
                status_text_main.text("QC..."); sc.pp.filter_cells(adata_pipeline, min_genes=st.session_state.min_genes_val); sc.pp.filter_genes(adata_pipeline, min_cells=st.session_state.min_cells_val)
                adata_pipeline.var['mt'] = adata_pipeline.var_names.str.startswith(st.session_state.mito_prefix_val); sc.pp.calculate_qc_metrics(adata_pipeline, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
                adata_pipeline = adata_pipeline[adata_pipeline.obs.pct_counts_mt < st.session_state.max_mito_pct_val, :].copy(); progress_bar_main.progress(12)
                
                status_text_main.text("Normalización..."); sc.pp.normalize_total(adata_pipeline, target_sum=1e4); sc.pp.log1p(adata_pipeline); progress_bar_main.progress(25)
                
                status_text_main.text("HVGs..."); sc.pp.highly_variable_genes(adata_pipeline, n_top_genes=st.session_state.n_top_genes_hvg_val, flavor='seurat_v3')
                if 'highly_variable' in adata_pipeline.var.columns and adata_pipeline.var.highly_variable.sum() > 0:
                    st.session_state.adata_hvg_filtered_intermediate = adata_pipeline[:, adata_pipeline.var.highly_variable].copy()
                else: raise ValueError("No se encontraron HVGs o la columna 'highly_variable' no se creó.")
                progress_bar_main.progress(37)

                if st.session_state.adata_hvg_filtered_intermediate is None: raise ValueError("adata_hvg_filtered_intermediate es None antes de Escalado.")
                
                status_text_main.text("Escalado..."); sc.pp.scale(st.session_state.adata_hvg_filtered_intermediate, max_value=10); progress_bar_main.progress(50)
                status_text_main.text("PCA..."); sc.tl.pca(st.session_state.adata_hvg_filtered_intermediate, svd_solver='arpack', n_comps=st.session_state.n_pcs_val); progress_bar_main.progress(62)
                status_text_main.text("Vecinos y UMAP..."); sc.pp.neighbors(st.session_state.adata_hvg_filtered_intermediate, n_neighbors=10, n_pcs=st.session_state.n_pcs_val); sc.tl.umap(st.session_state.adata_hvg_filtered_intermediate); progress_bar_main.progress(75)
                status_text_main.text("Clustering..."); 
                sc.tl.leiden(
                    st.session_state.adata_hvg_filtered_intermediate, 
                    resolution=st.session_state.leiden_res_val, 
                    key_added="leiden_clusters",
                    flavor="igraph",         # Usar el backend igraph
                    n_iterations=2,         # Como sugiere el warning para igraph
                    directed=False          # Asegurarse de que es False, como requiere igraph (generalmente es el default)
                )
                adata_pipeline.obs['leiden_clusters'] = st.session_state.adata_hvg_filtered_intermediate.obs['leiden_clusters'].copy()
                adata_pipeline.obsm['X_umap'] = st.session_state.adata_hvg_filtered_intermediate.obsm['X_umap'].copy()

                status_text_main.text("Marcadores de clúster..."); sc.tl.rank_genes_groups(adata_pipeline, 'leiden_clusters', method='wilcoxon', key_added='rank_genes_leiden_clusters')
                res_mk = adata_pipeline.uns['rank_genes_leiden_clusters']; grps_mk = res_mk['names'].dtype.names
                data_mk = [{'Cluster': grp, 'Rank': i_mk + 1, 'Gene': res_mk['names'][i_mk][grp], 'Score': res_mk['scores'][i_mk][grp], 'Log2FC': res_mk['logfoldchanges'][i_mk][grp], 'P-Value Adj': res_mk['pvals_adj'][i_mk][grp]} for grp in grps_mk for i_mk in range(st.session_state.n_top_markers_val) if i_mk < len(res_mk['names'][grp])]
                st.session_state.marker_genes_df = pd.DataFrame(data_mk); progress_bar_main.progress(100)
                
                status_text_main.empty(); st.session_state.adata_processed = adata_pipeline; st.session_state.analysis_done = True; st.balloons()
            except Exception as e_pipe: st.error(f"Error pipeline: {e_pipe}"); st.error(traceback.format_exc()); st.session_state.analysis_done = False; st.session_state.adata_processed = None
st.sidebar.markdown("---"); st.sidebar.info("App scRNA-seq v0.4")

# --- Sección de Resultados ---
st.header("Resultados del Análisis")
st.subheader("🔬 Explorador de Expresión Génica") # Input de genes siempre visible
st.session_state.gene_explorer_input = st.text_area("Ingresa nombres de genes (coma/espacio):", value=st.session_state.gene_explorer_input, key="gene_explorer_input_main_widget")

if st.session_state.analysis_done and st.session_state.adata_processed is not None:
    adata_res = st.session_state.adata_processed.copy()
    
    genes_to_viz_list = []
    genes_not_found_list = []
    if st.session_state.gene_explorer_input:
        raw_genes = [g.strip() for g in st.session_state.gene_explorer_input.replace(',', ' ').split()]
        genes_to_viz_list = [g for g in raw_genes if g and g in adata_res.var_names]
        genes_not_found_list = [g for g in raw_genes if g and g not in adata_res.var_names]
        if genes_not_found_list: st.warning(f"Genes no encontrados: {', '.join(genes_not_found_list)}")

    try: # Descarga AnnData
        with tempfile.NamedTemporaryFile(delete=False, suffix=".h5ad") as tmp_ad: path_tmp_ad = tmp_ad.name; adata_res.write_h5ad(path_tmp_ad)
        with open(path_tmp_ad, "rb") as f_ad: bytes_ad = f_ad.read()
        os.remove(path_tmp_ad)
        st.sidebar.download_button("Descargar AnnData Procesado (.h5ad)", bytes_ad, "processed_adata.h5ad", "application/octet-stream", key="download_adata_btn")
    except Exception as e_dl_ad: st.sidebar.error(f"Error descarga AnnData: {e_dl_ad}")

    tab_titles = ["📊 UMAPs", "🔬 Marcadores Clúster", "🧬 QC Plots", "📈 Análisis Diferencial", "🧬 Explorador Genes", "ℹ️ Info Dataset"]
    tab_umap, tab_markers, tab_qc, tab_dea, tab_gene_explorer, tab_info = st.tabs(tab_titles)

    with tab_umap:
        st.subheader("UMAP por Clústeres Leiden")
        fig_u_c, ax_u_c = plt.subplots(figsize=(7,6)); sc.pl.umap(adata_res,color='leiden_clusters',legend_loc='on data',ax=ax_u_c,show=False,title=f"Res: {st.session_state.leiden_res_val}"); st.pyplot(fig_u_c)
        st.download_button("UMAP Clúster (PNG)", fig_to_bytes(fig_u_c), "umap_clusters.png", "image/png", key="dl_u_c")
        plt.close(fig_u_c)

        st.subheader("UMAP por Muestra")
        if 'sample' in adata_res.obs:
            adata_res.obs['sample'] = adata_res.obs['sample'].astype('category')
            fig_u_s, ax_u_s = plt.subplots(figsize=(7,6)); sc.pl.umap(adata_res, color='sample', ax=ax_u_s, show=False, title="Por Muestra"); st.pyplot(fig_u_s)
            st.download_button("UMAP Muestra (PNG)", fig_to_bytes(fig_u_s), "umap_sample.png", "image/png", key="dl_u_s")
            plt.close(fig_u_s)
        
        st.subheader("UMAPs por Muestra (Coloreado por Clúster)")
        if 'sample' in adata_res.obs and 'leiden_clusters' in adata_res.obs:
            try:
                unique_s_list = sorted(adata_res.obs['sample'].cat.categories.tolist()); n_s = len(unique_s_list)
                cols_s_facet = min(n_s, 3); rows_s_facet = (n_s + cols_s_facet - 1) // cols_s_facet
                if n_s > 0:
                    fig_s_split, axes_s_split = plt.subplots(rows_s_facet, cols_s_facet, figsize=(cols_s_facet*5, rows_s_facet*4.5))
                    axes_s_flat = [axes_s_split] if n_s == 1 else axes_s_split.flatten()
                    for idx, s_val in enumerate(unique_s_list):
                        if idx < len(axes_s_flat):
                            ax_s_curr = axes_s_flat[idx]; tmp_ad_s = adata_res[adata_res.obs['sample'] == s_val].copy()
                            if not tmp_ad_s.obs.empty: sc.pl.umap(tmp_ad_s,color='leiden_clusters',ax=ax_s_curr,show=False,title=f"Muestra: {s_val}",legend_loc='on data' if idx==0 else None,legend_fontsize=6)
                            else: ax_s_curr.text(0.5,0.5,f"M: {s_val}\n(Sin células?)",ha='center',va='center'); ax_s_curr.set_xticks([]); ax_s_curr.set_yticks([])
                    for j_s in range(idx + 1, len(axes_s_flat)): fig_s_split.delaxes(axes_s_flat[j_s])
                    plt.tight_layout(); st.pyplot(fig_s_split)
                    st.download_button("UMAPs por Muestra (PNG)", fig_to_bytes(fig_s_split), "umaps_split_by_sample.png", "image/png", key="dl_u_split")
                    plt.close(fig_s_split)
            except Exception as e_u_split: st.error(f"Error UMAPs divididos: {e_u_split}"); st.error(traceback.format_exc())
    
    with tab_markers:
        st.subheader(f"Top {st.session_state.n_top_markers_val} Genes Marcadores por Clúster")
        if st.session_state.marker_genes_df is not None and not st.session_state.marker_genes_df.empty:
            st.dataframe(st.session_state.marker_genes_df)
            st.download_button("Marcadores Clúster (CSV)", st.session_state.marker_genes_df.to_csv(index=False).encode('utf-8'), "cluster_markers.csv", "text/csv", key="dl_markers_csv")
            
            st.subheader("Dot Plot Marcadores Clúster") # Lógica de dot plot para marcadores de clúster
            num_genes_dot_mk = min(5, st.session_state.n_top_markers_val); top_mk_genes_list = []
            if 'rank_genes_leiden_clusters' in adata_res.uns:
                for grp_mk_uns in adata_res.uns['rank_genes_leiden_clusters']['names'].dtype.names: top_mk_genes_list.extend(adata_res.uns['rank_genes_leiden_clusters']['names'][grp_mk_uns][:num_genes_dot_mk])
                top_mk_genes_plot = list(dict.fromkeys(top_mk_genes_list))
                if top_mk_genes_plot:
                    genes_in_ad_mk = [g for g in top_mk_genes_plot if g in adata_res.var_names]
                    if genes_in_ad_mk:
                        try:
                            fig_dot_mk, ax_dot_mk = plt.subplots(figsize=(max(10,len(genes_in_ad_mk)*0.4),max(5,len(adata_res.obs.leiden_clusters.cat.categories)*0.4)))
                            sc.pl.dotplot(adata_res,genes_in_ad_mk,groupby='leiden_clusters',ax=ax_dot_mk,show=False,standard_scale='var'); plt.xticks(rotation=90); plt.tight_layout(); st.pyplot(fig_dot_mk)
                            st.download_button("Dot Plot Marcadores (PNG)", fig_to_bytes(fig_dot_mk), "dotplot_cluster_markers.png", "image/png", key="dl_dot_mk")
                            plt.close(fig_dot_mk)
                        except Exception as e_dot_mk_plot: st.error(f"Error dot plot marcadores: {e_dot_mk_plot}"); st.error(traceback.format_exc())
                    else: st.warning("Genes para dot plot marcadores no encontrados.")
                else: st.info("No hay genes marcadores para dot plot.")
            else: st.info("Resultados de marcadores no disponibles para dot plot.")
        else: st.info("No se encontraron genes marcadores de clúster.")

    with tab_qc:
        st.subheader("QC Plots por Muestra") # Lógica de QC plots por muestra
        qc_metrics = ['n_genes_by_counts', 'total_counts', 'pct_counts_mt']; qc_titles_list = ["N Genes/Célula", "Total Cuentas/Célula", "% Cuentas Mito"]
        if 'sample' in adata_res.obs:
            adata_res.obs['sample'] = adata_res.obs['sample'].astype('category')
            for metric_qc, title_qc in zip(qc_metrics, qc_titles_list):
                st.markdown(f"#### {title_qc}")
                fig_qc_m, ax_qc_m = plt.subplots(figsize=(max(6,len(adata_res.obs['sample'].cat.categories)*1.5),5))
                try:
                    sc.pl.violin(adata_res,keys=metric_qc,groupby='sample',rotation=45,ax=ax_qc_m,show=False,cut=0); ax_qc_m.set_xlabel("Muestra"); plt.tight_layout(); st.pyplot(fig_qc_m)
                    st.download_button(f"Descargar {title_qc.replace('/','_')} (PNG)", fig_to_bytes(fig_qc_m), f"qc_violin_{metric_qc}_by_sample.png", "image/png", key=f"dl_qc_{metric_qc}")
                    plt.close(fig_qc_m)
                except Exception as e_qc_vln: st.error(f"Error violin QC para {metric_qc}: {e_qc_vln}"); plt.close(fig_qc_m)
        else: st.warning("Columna 'sample' no encontrada para QC plots comparativos.")

    with tab_dea:
        st.subheader("Análisis de Expresión Diferencial") # Lógica de DEA
        if st.session_state.dea_results_df is not None and not st.session_state.dea_results_df.empty:
            st.markdown(f"**Comparación:** `{st.session_state.dea_comparison_str}`")
            st.dataframe(st.session_state.dea_results_df.head(st.session_state.dea_n_genes_display))
            st.download_button("Tabla DEA Completa (CSV)", st.session_state.dea_results_df.to_csv(index=False).encode('utf-8'), f"dea_results_{st.session_state.dea_comparison_str.replace(' vs ','_vs_').replace(' ','_')}.csv", "text/csv", key="dl_dea_csv")
            
            st.markdown("#### Volcano Plot")
            try:
                df_volc = st.session_state.dea_results_df.copy(); df_volc['Significancia'] = 'No Significativo'
                df_volc.loc[(df_volc['Log2FC'] > st.session_state.dea_lfc_cutoff) & (df_volc['P-Value Adj'] < st.session_state.dea_pval_cutoff), 'Significancia'] = 'Upregulated'
                df_volc.loc[(df_volc['Log2FC'] < -st.session_state.dea_lfc_cutoff) & (df_volc['P-Value Adj'] < st.session_state.dea_pval_cutoff), 'Significancia'] = 'Downregulated'
                fig_volc_plot = px.scatter(df_volc,x='Log2FC',y='P-Value Adj',color='Significancia',color_discrete_map={'Upregulated':'red','Downregulated':'blue','No Significativo':'grey'},hover_data=['Gene','Log2FC','P-Value Adj'],title=f"Volcano: {st.session_state.dea_comparison_str}",labels={'Log2FC':'Log2 Fold Change','P-Value Adj':'P-valor Ajustado'})
                fig_volc_plot.update_layout(yaxis_type="log",yaxis_autorange="reversed")
                fig_volc_plot.add_hline(y=st.session_state.dea_pval_cutoff,line_dash="dash",line_color="black",annotation_text=f"P.adj={st.session_state.dea_pval_cutoff}")
                fig_volc_plot.add_vline(x=st.session_state.dea_lfc_cutoff,line_dash="dash",line_color="black",annotation_text=f"LFC={st.session_state.dea_lfc_cutoff}")
                fig_volc_plot.add_vline(x=-st.session_state.dea_lfc_cutoff,line_dash="dash",line_color="black",annotation_text=f"LFC={-st.session_state.dea_lfc_cutoff}")
                st.plotly_chart(fig_volc_plot, use_container_width=True)
                html_buf_volc = io.StringIO(); fig_volc_plot.write_html(html_buf_volc)
                st.download_button("Volcano Plot (HTML)", html_buf_volc.getvalue(), f"volcano_{st.session_state.dea_comparison_str.replace(' vs ','_vs_').replace(' ','_')}.html", "text/html", key="dl_volc_html")
            except Exception as e_volc_plot: st.error(f"Error Volcano Plot: {e_volc_plot}"); st.error(traceback.format_exc())
        elif st.session_state.get("dea_group1"): st.info("Ejecuta el Análisis Diferencial desde la barra lateral.")
        else: st.info("Configura el Análisis Diferencial en la barra lateral.")

    with tab_gene_explorer: # Lógica del explorador de genes
        st.subheader("Visualización de Genes Específicos")
        if not genes_to_viz_list: st.info("Ingresa nombres de genes válidos en el campo de arriba.")
        else:
            st.write(f"Visualizando: **{', '.join(genes_to_viz_list)}**")
            if genes_not_found_list: st.warning(f"No encontrados: {', '.join(genes_not_found_list)}")

            st.markdown("#### UMAPs por Expresión Génica") # UMAPs para genes del explorador
            n_g_u = len(genes_to_viz_list); c_g_u = min(n_g_u,3); r_g_u = (n_g_u+c_g_u-1)//c_g_u
            if n_g_u > 0:
                fig_ge_u, axes_ge_u = plt.subplots(r_g_u,c_g_u,figsize=(c_g_u*5,r_g_u*4.5))
                axes_ge_u_flat = [axes_ge_u] if n_g_u==1 else axes_ge_u.flatten()
                for i_ge_u,g_name_u in enumerate(genes_to_viz_list):
                    if i_ge_u < len(axes_ge_u_flat):
                        ax_g_u = axes_ge_u_flat[i_ge_u]
                        try: sc.pl.umap(adata_res,color=g_name_u,ax=ax_g_u,show=False,title=g_name_u,cmap='viridis',use_raw=False)
                        except: ax_g_u.text(0.5,0.5,f"Error\n{g_name_u}",ha='center',color='red'); ax_g_u.set_xticks([]); ax_g_u.set_yticks([])
                for j_ge_u in range(i_ge_u+1,len(axes_ge_u_flat)): fig_ge_u.delaxes(axes_ge_u_flat[j_ge_u])
                plt.tight_layout(); st.pyplot(fig_ge_u)
                st.download_button("UMAPs Genes (PNG)",fig_to_bytes(fig_ge_u),"gene_umaps.png","image/png",key="dl_ge_u")
                plt.close(fig_ge_u)

            st.markdown("#### Violines por Clúster") # Violines por clúster para genes del explorador
            if 'leiden_clusters' in adata_res.obs:
                n_g_v_cl = len(genes_to_viz_list); c_g_v_cl = min(n_g_v_cl,2); r_g_v_cl = (n_g_v_cl+c_g_v_cl-1)//c_g_v_cl
                if n_g_v_cl > 0:
                    fig_ge_v_cl,axes_ge_v_cl=plt.subplots(r_g_v_cl,c_g_v_cl,figsize=(c_g_v_cl*7,r_g_v_cl*5))
                    axes_ge_v_cl_flat = [axes_ge_v_cl] if n_g_v_cl==1 else axes_ge_v_cl.flatten()
                    for i_ge_v_cl,g_name_v_cl in enumerate(genes_to_viz_list):
                        if i_ge_v_cl < len(axes_ge_v_cl_flat):
                            ax_g_v_cl=axes_ge_v_cl_flat[i_ge_v_cl]
                            try: sc.pl.violin(adata_res,keys=g_name_v_cl,groupby='leiden_clusters',rotation=45,ax=ax_g_v_cl,show=False,use_raw=False,cut=0); ax_g_v_cl.set_title(g_name_v_cl)
                            except: ax_g_v_cl.text(0.5,0.5,f"Error\n{g_name_v_cl}",ha='center',color='red'); ax_g_v_cl.set_xticks([]); ax_g_v_cl.set_yticks([])
                    for j_ge_v_cl in range(i_ge_v_cl+1,len(axes_ge_v_cl_flat)): fig_ge_v_cl.delaxes(axes_ge_v_cl_flat[j_ge_v_cl])
                    plt.tight_layout(); st.pyplot(fig_ge_v_cl)
                    st.download_button("Violines por Clúster (PNG)",fig_to_bytes(fig_ge_v_cl),"gene_violins_cluster.png","image/png",key="dl_ge_v_cl")
                    plt.close(fig_ge_v_cl)
            
            if 'condition' in adata_res.obs and len(adata_res.obs['condition'].unique()) > 1: # Violines por condición para genes del explorador
                st.markdown("#### Violines por Condición")
                n_g_v_co = len(genes_to_viz_list); c_g_v_co = min(n_g_v_co,2); r_g_v_co = (n_g_v_co+c_g_v_co-1)//c_g_v_co
                if n_g_v_co > 0:
                    fig_ge_v_co,axes_ge_v_co=plt.subplots(r_g_v_co,c_g_v_co,figsize=(c_g_v_co*7,r_g_v_co*5))
                    axes_ge_v_co_flat = [axes_ge_v_co] if n_g_v_co==1 else axes_ge_v_co.flatten()
                    for i_ge_v_co,g_name_v_co in enumerate(genes_to_viz_list):
                        if i_ge_v_co < len(axes_ge_v_co_flat):
                            ax_g_v_co=axes_ge_v_co_flat[i_ge_v_co]
                            try: sc.pl.violin(adata_res,keys=g_name_v_co,groupby='condition',rotation=45,ax=ax_g_v_co,show=False,use_raw=False,cut=0); ax_g_v_co.set_title(g_name_v_co)
                            except: ax_g_v_co.text(0.5,0.5,f"Error\n{g_name_v_co}",ha='center',color='red'); ax_g_v_co.set_xticks([]); ax_g_v_co.set_yticks([])
                    for j_ge_v_co in range(i_ge_v_co+1,len(axes_ge_v_co_flat)): fig_ge_v_co.delaxes(axes_ge_v_co_flat[j_ge_v_co])
                    plt.tight_layout(); st.pyplot(fig_ge_v_co)
                    st.download_button("Violines por Condición (PNG)",fig_to_bytes(fig_ge_v_co),"gene_violins_condition.png","image/png",key="dl_ge_v_co")
                    plt.close(fig_ge_v_co)

            if len(genes_to_viz_list) > 1 and 'leiden_clusters' in adata_res.obs: # Dot plot para genes del explorador
                st.markdown("#### Dot Plot Genes Seleccionados por Clúster")
                try:
                    fig_ge_dot,ax_ge_dot=plt.subplots(figsize=(max(8,len(genes_to_viz_list)*0.7),max(5,len(adata_res.obs.leiden_clusters.cat.categories)*0.5)))
                    sc.pl.dotplot(adata_res,genes_to_viz_list,groupby='leiden_clusters',ax=ax_ge_dot,show=False,standard_scale='var',use_raw=False)
                    plt.xticks(rotation=90); plt.tight_layout(); st.pyplot(fig_ge_dot)
                    st.download_button("Dot Plot Genes (PNG)",fig_to_bytes(fig_ge_dot),"selected_genes_dotplot.png","image/png",key="dl_ge_dot")
                    plt.close(fig_ge_dot)
                except Exception as e_ge_dot_plot: st.error(f"Error dot plot genes: {e_ge_dot_plot}")

    with tab_info:
        st.subheader("Información del Dataset Procesado")
        st.write(f"Células: {adata_res.n_obs}, Genes: {adata_res.n_vars}")
        if 'sample' in adata_res.obs: st.write(f"Distribución por muestra:"); st.dataframe(adata_res.obs['sample'].value_counts())
        if st.session_state.adata_hvg_filtered_intermediate is not None: st.write(f"HVGs usados para PCA/etc.: {st.session_state.adata_hvg_filtered_intermediate.n_vars}")
        st.write("Metadatos células (obs) head:"); st.dataframe(adata_res.obs.head())
        st.write("Metadatos genes (var) head:"); st.dataframe(adata_res.var.head())
else:
    if st.session_state.adata_combined_raw is None and st.session_state.get("num_samples_value",0) > 0 and not all_files_uploaded_check :
        st.info("Sube los archivos para cada muestra y haz clic en 'Cargar y Concatenar Datos'.")
    elif st.session_state.get("num_samples_value",0) == 0 :
        st.info("Bienvenido. Ajusta el 'Número de muestras a cargar' para comenzar.")
    elif st.session_state.adata_combined_raw is None and st.session_state.get("num_samples_value",0) > 0 and all_files_uploaded_check:
        st.info("Archivos listos. Haz clic en 'Cargar y Concatenar Datos'.")
    elif st.session_state.adata_combined_raw is not None and not st.session_state.analysis_done:
        st.info("Datos cargados. Haz clic en 'Ejecutar Pipeline Principal'.")
    else: # Caso inicial sin acción tomada aún pero con num_samples > 0
        st.info("Bienvenido. Configura tus muestras y parámetros en la barra lateral.")
        

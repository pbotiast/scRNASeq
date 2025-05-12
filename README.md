# Analizador Interactivo de Single-Cell RNA-seq

## 1. Introducción

Este proyecto es una aplicación web interactiva desarrollada con Streamlit y Scanpy, diseñada para facilitar el análisis de datos de secuenciación de ARN de célula única (scRNA-seq). Permite a los usuarios cargar datos de múltiples muestras (en formato 10x Genomics), realizar un pipeline de análisis estándar que incluye control de calidad, normalización, reducción de dimensionalidad, clustering, y la identificación de genes marcadores. Adicionalmente, ofrece la posibilidad de realizar análisis de expresión diferencial entre condiciones definidas por el usuario.

**Propósito de la Aplicación:**
Proporcionar una herramienta visual e intuitiva para el análisis exploratorio de datos scRNA-seq, accesible tanto para biólogos con conocimientos básicos de bioinformática como para bioinformáticos que buscan una forma rápida de visualizar y procesar sus datos.

**Audiencia:**
Investigadores y científicos que trabajan con datos de scRNA-seq y necesitan una plataforma para realizar análisis preliminares y generar visualizaciones de forma interactiva.

**Tecnologías Clave:**
*   **Streamlit:** Framework para la creación de aplicaciones web interactivas con Python.
*   **Scanpy:** Biblioteca de Python para el análisis de datos de single-cell.
*   **AnnData:** Estructura de datos utilizada por Scanpy para almacenar matrices de expresión y metadatos.

## 2. Guía de Inicio Rápido

Sigue estos pasos para comenzar a analizar tus datos:

### Requisitos Previos:
*   Un navegador web moderno (Chrome, Firefox, Edge, Safari).
*   Tus archivos de datos de scRNA-seq en formato 10x Genomics. Para cada muestra, necesitarás:
    *   `matrix.mtx.gz` (o `matrix.mtx`)
    *   `features.tsv.gz` (o `genes.tsv.gz`, o `features.tsv`)
    *   `barcodes.tsv.gz` (o `barcodes.tsv`)

### Pasos para un Análisis Básico:

1.  **Ejecutar la Aplicación:**
    *   Si tienes el código localmente y Python/Streamlit instalados, abre una terminal o línea de comandos.
    *   Navega hasta el directorio donde guardaste el script (ej: `analizador_scRNAseq.py`).
    *   Ejecuta el comando: `streamlit run analizador_scRNAseq.py` (reemplaza `analizador_scRNAseq.py` con el nombre real de tu archivo).
    *   La aplicación se abrirá automáticamente en tu navegador web.

2.  **Configurar la Carga de Datos (en la barra lateral izquierda):**
    *   **Número de muestras a cargar:** Introduce cuántas muestras diferentes vas a analizar.
    *   Para cada muestra que aparezca:
        *   **Nombre Muestra X:** Escribe un nombre descriptivo para la muestra (ej: "Control_Dia0", "Tratado_Rep1").
        *   **Matrix.mtx (MX):** Haz clic en "Browse files" y selecciona el archivo `matrix.mtx.gz` (o `.mtx`) de esa muestra.
        *   **Features.tsv (MX):** Haz clic en "Browse files" y selecciona el archivo `features.tsv.gz` (o similar) de esa muestra.
        *   **Barcodes.tsv (MX):** Haz clic en "Browse files" y selecciona el archivo `barcodes.tsv.gz` (o `.tsv`) de esa muestra.

3.  **Cargar y Concatenar Datos:**
    *   Una vez que hayas subido todos los archivos para todas tus muestras, haz clic en el botón **"Cargar y Concatenar Datos"** en la barra lateral.
    *   Espera a que aparezca el mensaje de éxito (ej: "Cargado: X células, Y genes.").

4.  **Ajustar Parámetros del Pipeline (opcional):**
    *   En la sección "2. Parámetros de Pipeline Principal" de la barra lateral, puedes ajustar los valores para el filtrado, normalización, etc. Los valores por defecto son un buen punto de partida.

5.  **Ejecutar el Pipeline Principal:**
    *   Haz clic en el botón **"Ejecutar Pipeline Principal"** en la barra lateral.
    *   Este proceso puede tardar unos minutos dependiendo del tamaño de tus datos. Verás una barra de progreso y mensajes de estado.
    *   Una vez completado, aparecerá un mensaje de éxito (y unos globos 🎈).

6.  **Explorar los Resultados:**
    *   La sección principal de la aplicación se poblará con varias pestañas (`📊 UMAPs`, `🔬 Marcadores Clúster`, etc.). Navega por ellas para ver los resultados.
    *   Utiliza el **"🔬 Explorador de Expresión Génica"** en la parte superior de la sección de resultados para visualizar genes específicos.

7.  **Análisis de Expresión Diferencial (DEA) (opcional):**
    *   Si deseas comparar grupos de muestras, ve a la sección "3. Análisis de Expresión Diferencial (DEA)" en la barra lateral (esta sección aparece después de ejecutar el pipeline principal).
    *   **Asigna condiciones** a tus muestras (ej: Muestra1 -> "Control", Muestra2 -> "Tratado").
    *   **Selecciona Grupo 1 y Grupo 2** para la comparación.
    *   Ajusta los parámetros del DEA si es necesario.
    *   Haz clic en **"Ejecutar Análisis Diferencial"**.
    *   Los resultados aparecerán en la pestaña `📈 Análisis Diferencial`.

8.  **Descargar Resultados:**
    *   Utiliza los botones "Descargar..." que aparecen debajo de cada tabla o gráfico para guardar tus resultados.
    *   En la barra lateral, también encontrarás un botón para descargar el objeto `AnnData` procesado completo.

## 3. Descripción Detallada de la Interfaz y Funcionalidades

Esta sección describe en detalle cada componente de la interfaz de usuario y su función.

### 3.1. Barra Lateral (Sidebar)

La barra lateral, ubicada a la izquierda de la aplicación, contiene todos los controles para la carga de datos, la configuración de los parámetros del análisis y la ejecución de los principales pasos del pipeline.

#### 3.1.1. Sección "1. Carga de Datos"

Esta es la primera sección que encontrarás y es esencial para iniciar cualquier análisis.

*   **`Número de muestras a cargar`**:
    *   **Función:** Un campo numérico para especificar cuántas muestras biológicas o datasets individuales de scRNA-seq deseas analizar conjuntamente.
    *   **Uso:** Incrementa o disminuye este valor para que aparezcan los campos de carga de archivos correspondientes a cada muestra. Mínimo 1, máximo 10 (configurable en el código).

*   **Para cada Muestra (Muestra 1, Muestra 2, etc.):** (Organizado con `st.subheader` por muestra)
    *   **`Nombre Muestra X`**:
        *   **Función:** Un campo de texto para asignar un nombre identificador único a cada muestra. **Este nombre se conservará y se utilizará en todo el análisis** para identificar las células de esta muestra (ej. en plots UMAP, DEA, etc.).
        *   **Uso:** Introduce un nombre descriptivo (ej: `Control_Rep1`, `TratamientoA_Dia3`). Por defecto, se asigna `MuestraX`.
    *   **`Matrix (.mtx/.mtx.gz)`**:
        *   **Función:** Botón para subir el archivo de matriz de cuentas.
        *   **Uso:** Selecciona el archivo `.mtx` o `.mtx.gz` correspondiente.
    *   **`Features (.tsv/.tsv.gz)`**:
        *   **Función:** Botón para subir el archivo de características/genes.
        *   **Uso:** Selecciona el archivo `.tsv` o `.tsv.gz` (o `genes.tsv`/`genes.tsv.gz`) correspondiente.
    *   **`Barcodes (.tsv/.tsv.gz)`**:
        *   **Función:** Botón para subir el archivo de códigos de barras celulares.
        *   **Uso:** Selecciona el archivo `.tsv` o `.tsv.gz` correspondiente.

*   **Botón `Cargar y Concatenar Datos`**:
    *   **Función:**
        1.  **Validación de Archivos:** Antes de cargar, la aplicación realiza una validación básica del formato de los archivos subidos (ej. si el matrix.mtx parece un archivo MatrixMarket). Los resultados de la validación se muestran en un expander.
        2.  **Carga y Concatenación:** Si todos los archivos son válidos, inicia el proceso de carga. Los datos de cada muestra se leen individualmente (conservando el nombre de muestra proporcionado) y luego se concatenan en un único objeto AnnData (`adata_raw`). La columna `adata_raw.obs['sample']` contendrá los nombres de muestra que especificaste.
    *   **Uso:** Púlsalo *después* de haber seleccionado todos los archivos necesarios para todas las muestras.
    *   **Nota:** Si la validación falla para alguna muestra, la carga no procederá. Deberás corregir los archivos y volver a intentarlo. Al pulsar este botón, se reiniciarán los resultados de cualquier pipeline anterior.

#### 3.1.2. Sección "2. Parámetros del Pipeline"

Esta sección te permite configurar los parámetros para los pasos de preprocesamiento, análisis y clustering. El orden principal de las operaciones del pipeline es: QC -> HVG -> Normalización/Log -> Escalado (de HVGs) -> PCA -> Vecinos -> UMAP -> Leiden.

*   **`Mínimo genes/célula`**: Filtra células con un número de genes detectados inferior a este umbral.
*   **`Mínimo células/gen`**: Filtra genes que se expresan en un número de células inferior a este umbral.
*   **`Prefijo genes mitocondriales`**: Cadena para identificar genes mitocondriales (ej: `MT-` para humano).
*   **`Máx % cuentas mitocondriales`**: Porcentaje máximo de cuentas mitocondriales permitido por célula.
*   **`Nº HVGs a seleccionar`**: Número de Genes Altamente Variables (HVGs) a seleccionar. La selección de HVGs (usando el método `seurat_v3` con `batch_key='sample'`) se realiza **antes** de la normalización global, sobre los datos de cuentas crudas post-QC.
*   **`Nº PCs (para PCA y Vecinos)`**: Número de componentes principales a calcular con PCA y a usar para la construcción del grafo de vecinos. El valor se ajusta automáticamente si es demasiado alto para las dimensiones de los datos post-HVG. Mínimo 5 recomendado para el solver `arpack`.
*   **`Nº Vecinos (para grafo KNN)`**: Número de vecinos a considerar al construir el grafo de Vecinos Próximos (KNN), que se usa para UMAP y Leiden. Se ajusta automáticamente si es demasiado alto para el número de células.
*   **`Resolución Leiden`**: Parámetro del algoritmo de clustering Leiden. Valores más altos tienden a producir más clústeres.
*   **`Nº marcadores a mostrar/clúster`**: Cuántos genes marcadores se mostrarán en la tabla de resultados.
*   **`Backend Leiden`**: Permite elegir el backend para el algoritmo de Leiden (`igraph` o `leidenalg`). `igraph` es generalmente recomendado y es el default.

*   **Subsección `Parámetros UMAP`**:
    *   **`Calcular también UMAP 3D`**: Checkbox para opcionalmente calcular y permitir la visualización de un embedding UMAP en 3 dimensiones.
    *   **`Inicialización UMAP`**: Método de inicialización para UMAP (`spectral`, `random`, `pca`). `'random'` es el default actual en la app por mayor estabilidad con algunas combinaciones de versiones de bibliotecas, aunque `'spectral'` es a menudo preferido.
    *   **`Nº Vecinos UMAP (para embedding)`**: Número de vecinos que UMAP considera al construir su propia representación del grafo para el embedding. Este parámetro es específico de UMAP y puede ser diferente del "Nº Vecinos (para grafo KNN)". Controla el balance entre detalle local y estructura global en el plot UMAP.
    *   **`Distancia Mínima UMAP`**: Controla cuán agrupados o dispersos estarán los puntos en el embedding UMAP. Valores más bajos tienden a agrupar más los puntos.

*   **Subsección (Nueva) `Personalización de Plots`** (ubicada al final de la sidebar o en su propio expander):
    *   **`Paleta de Colores (Clusters/Muestras)`**: Permite seleccionar una paleta de colores de Matplotlib/Scanpy para los plots UMAP y otros.
    *   **`Tamaño de Puntos UMAP (aprox.)`**: Controla el tamaño de los puntos en los plots UMAP 2D generados con `sc.pl.umap`.
    *   **`Nº genes/clúster para Heatmap Marcadores`**: Define cuántos de los top marcadores por clúster se incluirán en la visualización del heatmap.

*   **Botón `Ejecutar Pipeline Principal`**:
    *   Inicia la secuencia completa de análisis con los parámetros configurados. Los pasos principales son:
        1.  Control de Calidad (QC).
        2.  Selección de HVGs (sobre datos crudos post-QC, usando `batch_key`).
        3.  Normalización y transformación logarítmica del dataset completo.
        4.  Creación de un subconjunto de datos conteniendo solo los HVGs (ya normalizados/log).
        5.  Escalado de este subconjunto de HVGs.
        6.  PCA y cálculo de Vecinos KNN (sobre el subconjunto de HVGs escalado).
        7.  Cálculo de UMAP 2D (y opcionalmente 3D) usando la API directa de `umap-learn` sobre los PCs.
        8.  Clustering con Leiden sobre el grafo KNN.
        9.  Transferencia de resultados (clusters, UMAPs) al AnnData procesado completo.
        10. Cálculo de genes marcadores.
    *   Si tiene éxito, las pestañas de resultados se actualizan.

#### 3.1.3. Sección "3. Análisis de Expresión Diferencial (DEA)"

Esta sección aparece en la barra lateral **únicamente después de que el "Pipeline Principal" se haya ejecutado con éxito** (es decir, `st.session_state.analysis_done` es `True` y `st.session_state.adata_processed` existe). Permite comparar la expresión génica entre diferentes grupos de muestras, que se definen como "condiciones".

*   **`Asignar Muestras a Condiciones`**:
    *   **Función:** Para cada muestra cargada originalmente (identificada por su nombre), se muestra un campo de texto donde puedes asignar una "condición" o grupo experimental.
    *   **Uso:** Escribe el nombre de la condición para cada muestra. Por ejemplo, si tienes `MuestraA` y `MuestraB` que son controles, y `MuestraC` y `MuestraD` que son tratadas, podrías asignar:
        *   `MuestraA`: "Control"
        *   `MuestraB`: "Control"
        *   `MuestraC`: "Tratado"
        *   `MuestraD`: "Tratado"
    *   **Importancia:** Debes definir al menos dos condiciones diferentes para poder realizar una comparación. Los nombres de las condiciones son sensibles a mayúsculas y minúsculas.

*   **`Seleccionar Grupos para Comparación`**:
    *   Esta subsección aparece si has definido al menos dos condiciones válidas.
    *   **`Grupo 1 (Referencia)`**:
        *   **Función:** Un menú desplegable para seleccionar la condición que servirá como grupo de referencia o control en la comparación.
        *   **Uso:** Elige una de las condiciones únicas que definiste arriba.
    *   **`Grupo 2 (Comparación)`**:
        *   **Función:** Un menú desplegable para seleccionar la condición que se comparará con el Grupo 1. Los genes sobreexpresados en el Grupo 2 (relativo al Grupo 1) tendrán un Log2FC positivo.
        *   **Uso:** Elige una condición diferente al Grupo 1.

*   **`Realizar DEA en:`**:
    *   **Función:** Un menú desplegable para especificar el alcance del análisis diferencial:
        *   **`Todos los Clústeres`**: El DEA se realizará utilizando todas las células pertenecientes a las muestras asignadas al Grupo 1 y Grupo 2, independientemente de su clúster de Leiden.
        *   **`Clúster X` (ej: Clúster 0, Clúster 1, etc.)**: El DEA se restringirá a las células que pertenecen a este clúster específico *y* a las muestras asignadas al Grupo 1 y Grupo 2. Esto es útil para encontrar diferencias específicas de un tipo celular entre condiciones.
    *   **Uso:** Selecciona el alcance deseado.

*   **`Nº genes DEA`**:
    *   **Función:** Un control deslizante para determinar cuántos de los genes más significativos (ordenados por P-valor ajustado) se mostrarán en la tabla de resultados en la pestaña "📈 Análisis Diferencial".
    *   **Rango:** 10-200 (configurable).
    *   **Nota:** La descarga CSV siempre contendrá todos los genes.

*   **`Log2FC cutoff`**:
    *   **Función:** Un campo numérico para establecer el umbral mínimo (en valor absoluto) del Log2 Fold Change (cambio logarítmico de la expresión) para que un gen sea considerado biológicamente significativo en el Volcano Plot. No filtra la tabla directamente, pero se usa para colorear el Volcano Plot.
    *   **Uso:** Introduce un valor (ej: 0.25, 0.5, 1.0). Un valor de 0.5 significa que se resaltarán genes que cambian al menos ~1.4 veces.

*   **`P-adj cutoff`**:
    *   **Función:** Un campo numérico para establecer el umbral del p-valor ajustado (corrected p-value) para que un gen sea considerado estadísticamente significativo en el Volcano Plot. No filtra la tabla directamente, pero se usa para colorear el Volcano Plot.
    *   **Uso:** Introduce un valor (ej: 0.05, 0.01).

*   **Botón `Ejecutar Análisis Diferencial`**:
    *   **Función:** Inicia el cálculo del DEA basado en los parámetros seleccionados.
        1.  Se crea una copia del AnnData procesado.
        2.  Se añade una columna `condition` a los metadatos de las células (`.obs`) basada en las asignaciones.
        3.  Se filtran las células para incluir solo aquellas pertenecientes al Grupo 1 y Grupo 2.
        4.  Si se seleccionó un clúster específico, se filtran adicionalmente las células para ese clúster.
        5.  Se ejecuta `scanpy.tl.rank_genes_groups` usando el método Wilcoxon para comparar el Grupo 2 contra el Grupo 1.
        6.  Los resultados se almacenan y se formatean en un DataFrame.
    *   **Uso:** Púlsalo después de configurar todas las opciones del DEA.
    *   **Resultado:** Si tiene éxito, los resultados se mostrarán en la pestaña "📈 Análisis Diferencial". Si hay muy pocas células en alguno de los grupos (menos de 3), se mostrará un error.

#### 3.1.4. Botones de Acción y Descarga Adicionales en la Sidebar

Al final de la barra lateral, o distribuidos en ella, pueden aparecer otros elementos:

*   **Mensajes de estado o advertencias:** Por ejemplo, un aviso si no se han subido todos los archivos para las muestras seleccionadas.
*   **`Descargar AnnData Procesado (.h5ad)`:** (Generalmente aparece después de ejecutar el pipeline principal).
    *   **Función:** Permite descargar el objeto AnnData completo (`adata_processed`) que contiene los datos crudos (si se mantuvieron), los datos normalizados, los metadatos de células y genes, los resultados del PCA, las coordenadas UMAP, las asignaciones de clústeres, y los resultados de los genes marcadores de clúster.
    *   **Uso:** Útil para análisis posteriores fuera de la aplicación, usando Scanpy en un entorno Python, o para compartir los datos procesados. El archivo está en formato HDF5.
*   **Información de la App:**
    *   Al final de la sidebar, se muestra la versión de la aplicación (ej: `App scRNA-seq v0.4`).

#### 3.1.5. Sección "Guardar/Cargar Configuración"

Esta sección permite guardar y cargar los parámetros de configuración del pipeline para facilitar la reproducibilidad y la aplicación de configuraciones consistentes a diferentes análisis.

*   **Botón `Guardar Configuración Actual`**:
    *   Al pulsarlo, se genera un archivo JSON que contiene los valores actuales de todos los parámetros configurables en la sidebar (excepto los datos AnnData en sí mismos y los archivos subidos).
    *   Se ofrece un botón para descargar este archivo `scRNAseq_app_params_[fecha].json`.
*   **`Cargar Configuración (.json)`**:
    *   Permite subir un archivo JSON previamente guardado.
    *   Si el archivo es válido, los parámetros en la sidebar se actualizarán con los valores del archivo.
    *   Es útil para restaurar una configuración de análisis anterior o para compartir parámetros.

### 3.2. Sección de Resultados (Panel Principal)

Una vez que el "Pipeline Principal" ha sido ejecutado, el panel principal de la aplicación, a la derecha de la barra lateral, se activa y muestra los resultados organizados en varias pestañas. En la parte superior de esta sección, siempre visible, se encuentra el "Explorador de Expresión Génica".

#### 3.2.1. 🔬 Explorador de Expresión Génica

Este campo de texto te permite visualizar rápidamente la expresión de genes específicos de tu interés a través de varios tipos de gráficos en la pestaña "🧬 Explorador Genes".

*   **`Ingresa nombres de genes (coma/espacio):`**:
    *   **Función:** Un área de texto donde puedes escribir o pegar una lista de nombres de genes. Los nombres de los genes deben coincidir con los presentes en tus datos (generalmente símbolos de genes).
    *   **Uso:** Separa los nombres de los genes con comas (`,`) o espacios. Por ejemplo: `CD4, CD8A, MS4A1, GAPDH`.
    *   **Feedback:**
        *   Si introduces genes que no se encuentran en el dataset, aparecerá una advertencia listándolos.
        *   Los gráficos en la pestaña "🧬 Explorador Genes" se generarán solo para los genes válidos encontrados.

#### 3.2.2. Pestaña: "📊 UMAPs"

Esta pestaña se centra en las visualizaciones UMAP (Uniform Manifold Approximation and Projection), que son proyecciones 2D de tus datos que intentan preservar la estructura global y local de las poblaciones celulares.

*   **`UMAP por Clústeres Leiden`**:
    *   **Visualización:** Un gráfico UMAP donde cada punto representa una célula, y el color de cada punto indica el clúster de Leiden al que ha sido asignada.
    *   **Interpretación:** Ayuda a visualizar la separación de los diferentes clústeres celulares identificados. Idealmente, células del mismo clúster deberían agruparse en el UMAP. Se incluye la resolución de Leiden utilizada.
    *   **Descarga:** Botón "UMAP Clúster (PNG)" para descargar la imagen.

*   **`UMAP por Muestra`**:
    *   **Visualización:** Un gráfico UMAP similar al anterior, pero esta vez el color de cada punto indica la muestra de origen de la célula (según los nombres que proporcionaste durante la carga de datos).
    *   **Interpretación:** Útil para identificar si hay efectos de batch (si las células de una muestra se agrupan separadamente de otras sin una razón biológica clara) o si ciertos tipos celulares son específicos de una muestra/condición.
    *   **Descarga:** Botón "UMAP Muestra (PNG)" para descargar la imagen.

*   **`UMAPs por Muestra (Coloreado por Clúster)`**:
    *   **Visualización:** Una serie de gráficos UMAP, uno por cada muestra cargada. Dentro de cada UMAP individual (correspondiente a una muestra), las células están coloreadas según su clúster de Leiden.
    *   **Interpretación:** Permite ver la distribución de los clústeres celulares dentro de cada muestra individualmente. Es muy útil para comparar cómo se componen las poblaciones celulares entre diferentes muestras o condiciones. Por ejemplo, podrías observar si un clúster específico está enriquecido o ausente en una muestra particular.
    *   **Formato:** Los gráficos se organizan en una cuadrícula para facilitar la comparación.
    *   **Descarga:** Botón "UMAPs por Muestra (PNG)" para descargar la imagen combinada.

#### 3.2.3. Pestaña: "🔬 Marcadores Clúster"

Esta pestaña muestra información sobre los genes que son característicos de cada clúster de Leiden identificado. Estos genes marcadores ayudan a inferir la identidad o el estado de los tipos celulares representados por cada clúster.

*   **`Top X Genes Marcadores por Clúster` (Tabla)**:
    *   **Visualización:** Una tabla que lista los `X` genes más significativos (según el parámetro `Nº marcadores/clúster` de la sidebar) para cada clúster de Leiden.
    *   **Columnas:**
        *   `Cluster`: El identificador del clúster de Leiden.
        *   `Rank`: La posición del gen dentro de los marcadores de ese clúster (ej: 1º, 2º, etc.).
        *   `Gene`: El nombre del gen marcador.
        *   `Score`: Puntuación del test estadístico (Wilcoxon) utilizado para identificar el marcador. Valores más altos indican una mayor especificidad.
        *   `Log2FC`: Log2 Fold Change. Indica cuánto más (o menos) se expresa el gen en el clúster actual en comparación con el promedio del resto de clústeres.
        *   `P-Value Adj`: P-valor ajustado por múltiples comparaciones. Indica la significancia estadística del gen como marcador.
    *   **Interpretación:** Los genes con Log2FC alto y P-valor ajustado bajo son buenos candidatos a marcadores específicos del clúster.
    *   **Descarga:** Botón "Marcadores Clúster (CSV)" para descargar la tabla completa en formato CSV.

*   **`Dot Plot Marcadores Clúster`**:
    *   **Visualización:** Un gráfico de puntos (dot plot) que muestra la expresión de un subconjunto de los genes marcadores principales (generalmente los 1-5 mejores de cada clúster, dependiendo del parámetro `Nº marcadores/clúster`) a través de todos los clústeres de Leiden.
    *   **Interpretación:**
        *   **Color del punto:** Indica el nivel promedio de expresión del gen en las células de ese clúster (generalmente, colores más cálidos o intensos significan mayor expresión).
        *   **Tamaño del punto:** Indica el porcentaje de células dentro de ese clúster que expresan el gen (puntos más grandes significan que más células del clúster expresan ese gen).
    *   Este gráfico es muy útil para visualizar patrones de expresión y confirmar la especificidad de los marcadores.
    *   **Descarga:** Botón "Dot Plot Marcadores (PNG)" para descargar la imagen.

#### 3.2.4. Pestaña: "🔥 Heatmap Marcadores"

Esta pestaña visualiza la expresión de los genes marcadores más importantes a través de los clústeres de Leiden en forma de heatmap.

*   **`Heatmap de Top X Genes Marcadores por Clúster`**:
    *   **Visualización:** Un heatmap donde las filas suelen ser genes y las columnas son células (agrupadas y promediadas por clúster, o mostrando células individuales). El color indica el nivel de expresión.
    *   **Selección de Genes:** Utiliza los "Nº genes/clúster para Heatmap Marcadores" definidos en la sidebar para seleccionar los N mejores marcadores de cada clúster.
    *   **Dendrograma:** Si se calcula con éxito (basándose en `X_pca_hvg`), se muestra un dendrograma que agrupa los clústeres según la similitud de su perfil de expresión para los genes mostrados.
    *   **Escalado:** La expresión suele estar escalada por gen (Z-score) para resaltar patrones relativos.
    *   **Interpretación:** Ayuda a ver patrones de co-expresión y la especificidad de los marcadores de forma visual.
    *   **Descarga:** Botón para descargar el heatmap como imagen PNG.
    *   
#### 3.2.5. Pestaña: "🧬 QC Plots"

Esta pestaña muestra gráficos de control de calidad (Quality Control) que resumen métricas importantes sobre las células, agrupadas por la muestra original. Estos gráficos se generan sobre los datos *después* del filtrado inicial.

*   **`QC Plots por Muestra`**:
    *   **Visualización:** Se generan varios diagramas de violín, uno para cada métrica de QC principal:
        *   **`N Genes/Célula`**: Distribución del número de genes detectados por célula, para cada muestra.
        *   **`Total Cuentas/Célula`**: Distribución del número total de transcritos (UMIs) contados por célula, para cada muestra.
        *   **`% Cuentas Mito`**: Distribución del porcentaje de cuentas provenientes de genes mitocondriales por célula, para cada muestra.
    *   **Interpretación:** Estos gráficos permiten comparar la calidad de las células entre diferentes muestras. Diferencias grandes podrían indicar problemas técnicos en alguna muestra o diferencias biológicas intrínsecas. Por ejemplo, después del filtrado, se espera que los porcentajes mitocondriales sean bajos y relativamente homogéneos.
    *   **Descarga:** Cada gráfico de violín tiene su propio botón de descarga "Descargar [NombreMétrica] (PNG)".

#### 3.2.6. Pestaña: "📈 Análisis Diferencial"

Esta pestaña muestra los resultados del Análisis de Expresión Diferencial (DEA) si se ha ejecutado desde la barra lateral.

*   **Información de la Comparación**:
    *   Se muestra un texto indicando qué grupos se están comparando (ej: `Comparación: Tratado vs Control (Clúster 3)`).

*   **Tabla de Resultados DEA**:
    *   **Visualización:** Una tabla que muestra los genes diferencialmente expresados entre los dos grupos seleccionados, ordenados por P-valor ajustado. El número de genes mostrados está determinado por el parámetro `Nº genes DEA` de la sidebar.
    *   **Columnas:**
        *   `Gene`: Nombre del gen.
        *   `Log2FC`: Log2 Fold Change. Un valor positivo indica que el gen está sobreexpresado en el "Grupo 2 (Comparación)" relativo al "Grupo 1 (Referencia)". Un valor negativo indica subexpresión.
        *   `P-Value`: P-valor crudo del test estadístico.
        *   `P-Value Adj`: P-valor ajustado por múltiples comparaciones (ej: Benjamini-Hochberg). Es el valor principal para determinar significancia.
        *   `Scores`: Puntuación del test (ej: estadístico U de Wilcoxon).
    *   **Descarga:** Botón "Tabla DEA Completa (CSV)" para descargar la tabla completa con todos los genes analizados.

*   **`Volcano Plot`**:
    *   **Visualización:** Un gráfico de dispersión interactivo (Plotly) que visualiza los resultados del DEA.
        *   **Eje X:** Log2 Fold Change.
        *   **Eje Y:** -log10(P-valor Ajustado). Valores más altos en este eje indican mayor significancia estadística.
    *   **Interpretación:**
        *   Cada punto es un gen.
        *   Los genes se colorean según su significancia:
            *   `Upregulated` (Sobreexpresados): Generalmente en rojo, son los genes con Log2FC positivo por encima del umbral y P-valor ajustado por debajo del umbral.
            *   `Downregulated` (Subexpresados): Generalmente en azul, son los genes con Log2FC negativo por debajo del umbral (más negativo que -umbral) y P-valor ajustado por debajo del umbral.
            *   `No Significativo`: Generalmente en gris.
        *   Se dibujan líneas que indican los umbrales de `Log2FC cutoff` y `P-adj cutoff` seleccionados en la sidebar.
        *   Al pasar el ratón sobre un punto, se muestra información adicional del gen.
    *   **Descarga:** Botón "Volcano Plot (HTML)" para descargar el gráfico interactivo como un archivo HTML independiente.

#### 3.2.7. Pestaña: "🧬 Explorador Genes"

Esta pestaña muestra visualizaciones específicas para los genes que has introducido en el campo "🔬 Explorador de Expresión Génica" en la parte superior del panel de resultados.

*   **Mensajes de Estado**:
    *   Si no se han introducido genes, se mostrará un mensaje indicándolo.
    *   Se listarán los genes que se están visualizando y aquellos que no se encontraron en el dataset.

*   **`UMAPs por Expresión Génica`**:
    *   **Visualización:** Una serie de gráficos UMAP, uno por cada gen válido introducido. En cada UMAP, las células se colorean según el nivel de expresión del gen correspondiente (generalmente usando un gradiente de color donde la intensidad indica mayor expresión).
    *   **Interpretación:** Permite ver en qué poblaciones celulares (según la distribución UMAP) se expresa cada gen de interés.
    *   **Descarga:** Botón "UMAPs Genes (PNG)" para descargar la imagen combinada.

*   **`Violines por Clúster`**:
    *   **Visualización:** Para cada gen válido introducido, se muestra un diagrama de violín que ilustra la distribución de su expresión a través de los diferentes clústeres de Leiden.
    *   **Interpretación:** Ayuda a cuantificar y comparar los niveles de expresión de un gen entre los clústeres identificados.
    *   **Descarga:** Botón "Violines por Clúster (PNG)" para descargar la imagen combinada.

*   **`Violines por Condición` (si las condiciones están definidas en el DEA)**:
    *   **Visualización:** Si se ha realizado un DEA y se han definido condiciones, para cada gen válido introducido se muestra un diagrama de violín que ilustra la distribución de su expresión a través de las diferentes condiciones asignadas a las muestras.
    *   **Interpretación:** Permite comparar los niveles de expresión de un gen entre las condiciones experimentales.
    *   **Descarga:** Botón "Violines por Condión (PNG)" para descargar la imagen combinada.

*   **`Dot Plot Genes Seleccionados por Clúster` (si se introduce más de un gen)**:
    *   **Visualización:** Similar al "Dot Plot Marcadores Clúster", pero utilizando los genes que el usuario ha introducido en el explorador.
    *   **Interpretación:** Muestra el nivel promedio de expresión (color) y el porcentaje de células que expresan (tamaño del punto) cada uno de los genes seleccionados, a través de todos los clústeres de Leiden.
    *   **Descarga:** Botón "Dot Plot Genes (PNG)" para descargar la imagen.

#### 3.2.7. Pestaña: "ℹ️ Info Dataset"

Esta pestaña proporciona información resumida y metadatos sobre el conjunto de datos procesado.

*   **Estadísticas Generales**:
    *   `Células`: Número total de células después del filtrado.
    *   `Genes`: Número total de genes después del filtrado.
    *   `Distribución por muestra`: Tabla o serie que muestra cuántas células provienen de cada muestra original.
    *   `HVGs usados para PCA/etc.`: Número de genes altamente variables que se utilizaron para los pasos de reducción de dimensionalidad y clustering.
*   **Vistazo a los Metadatos**:
    *   `Metadatos células (obs) head`: Muestra las primeras filas de la tabla de metadatos de las células (`adata.obs`). Esto incluye columnas como `sample`, `n_genes_by_counts`, `total_counts`, `pct_counts_mt`, `leiden_clusters`, y `condition` (si se hizo DEA).
    *   `Metadatos genes (var) head`: Muestra las primeras filas de la tabla de metadatos de los genes (`adata.var`). Esto puede incluir columnas como `n_cells_by_counts`, `mt` (marcador mitocondrial), `highly_variable`.

## 4. Formato de Datos de Entrada

Para utilizar el Analizador Interactivo de Single-Cell RNA-seq, necesitas proporcionar los datos de salida de un pipeline de preprocesamiento estándar de 10x Genomics (como Cell Ranger) para cada muestra. La aplicación espera tres archivos específicos por muestra, que generalmente se encuentran en el subdirectorio `filtered_feature_bc_matrix` (o similar) de la salida de Cell Ranger:

1.  **`matrix.mtx.gz` (o `matrix.mtx`) - Archivo de Matriz:**
    *   **Descripción:** Este archivo contiene la matriz de cuentas de expresión en formato Matrix Market. Es un formato de texto disperso que lista las cuentas no nulas.
    *   **Contenido:**
        *   La primera línea suele ser un comentario que empieza con `%%MatrixMarket`.
        *   La segunda línea (a veces también comentario) puede indicar el tipo de datos.
        *   La tercera línea contiene tres números: número de genes, número de códigos de barras (células), y número total de entradas no nulas en la matriz.
        *   Las líneas subsiguientes contienen tres valores por línea: `índice_gen índice_código_de_barras cuenta`. Los índices están basados en 1.
    *   **Compresión:** La aplicación puede manejar tanto la versión comprimida (`.gz`) como la no comprimida.

2.  **`features.tsv.gz` (o `genes.tsv.gz`, o `features.tsv`) - Archivo de Características/Genes:**
    *   **Descripción:** Este archivo de texto tabulado lista los genes (u otras características) cuyas cuentas se encuentran en el archivo `matrix.mtx`. El orden de los genes en este archivo corresponde al primer índice (índice de gen) en el archivo `matrix.mtx`.
    *   **Contenido Típico (puede variar ligeramente según la versión de Cell Ranger):**
        *   Generalmente dos columnas, separadas por tabulador:
            1.  ID del gen (ej: Ensembl ID como `ENSG00000243485`).
            2.  Símbolo del gen (ej: `MIR1302-2HG`). A veces, puede haber una tercera columna con el tipo de característica (ej: "Gene Expression").
        *   La aplicación utiliza la columna de símbolos de gen (o la que `scanpy.read_10x_mtx` interprete como `gene_symbols` al usar `var_names='gene_symbols'`).
    *   **Compresión:** Puede estar comprimido (`.gz`) o no.

3.  **`barcodes.tsv.gz` (o `barcodes.tsv`) - Archivo de Códigos de Barras:**
    *   **Descripción:** Este archivo de texto tabulado lista los códigos de barras de las células. El orden de los códigos de barras en este archivo corresponde al segundo índice (índice de código de barras) en el archivo `matrix.mtx`.
    *   **Contenido:**
        *   Una única columna con los códigos de barras de las células (ej: `AAACCCAAGGAGAGTA-1`).
    *   **Compresión:** Puede estar comprimido (`.gz`) o no.

**Importante:**
*   Asegúrate de que para cada muestra los tres archivos (`matrix`, `features`/`genes`, `barcodes`) provengan del mismo análisis y sean consistentes entre sí.
*   Los nombres de los archivos dentro del widget de carga deben corresponder a los nombres esperados por `scanpy.read_10x_mtx` cuando se le pasa un directorio (es decir, `matrix.mtx.gz`, `features.tsv.gz`, `barcodes.tsv.gz`). Aunque la aplicación maneja la carga de archivos individuales, el proceso interno los coloca en un directorio temporal para que Scanpy los lea. Si tus archivos tienen nombres ligeramente diferentes (ej: `genes.tsv` en lugar de `features.tsv`), la función `sc.read_10x_mtx` es generalmente robusta para encontrarlos, pero usar los nombres estándar es la práctica más segura.

## 5. Solución de Problemas (FAQ)

Aquí encontrarás respuestas a preguntas frecuentes y soluciones a problemas comunes que podrías encontrar al usar la aplicación.

**P1: La aplicación muestra un error al intentar "Cargar y Concatenar Datos". ¿Qué puede estar pasando?**
*   **Archivos incorrectos o corruptos:** Asegúrate de que has subido los archivos correctos (`matrix.mtx.gz`, `features.tsv.gz`, `barcodes.tsv.gz` o sus variantes no comprimidas) para cada muestra y que no están corruptos. Intenta abrirlos o descomprimirlos localmente para verificar.
*   **Formato inesperado dentro de los archivos:** Aunque `scanpy.read_10x_mtx` es robusto, desviaciones significativas del formato estándar 10x pueden causar problemas. Verifica que la estructura interna de tus archivos (especialmente las primeras líneas de `matrix.mtx` y el número de columnas en `features.tsv`) sea la esperada.
*   **Nombres de archivo no estándar:** Si bien la aplicación permite subir archivos con cualquier nombre, internamente los renombra o los coloca en una estructura donde Scanpy espera los nombres estándar (`matrix.mtx`, `features.tsv`, `barcodes.tsv`). Si tus archivos originales dentro de un posible directorio comprimido (que no es el caso aquí, ya que se suben individualmente) fueran muy diferentes, podría ser un problema para algunas funciones de lectura menos flexibles, pero la implementación actual debería ser robusta para archivos individuales.
*   **Inconsistencia entre archivos:** El número de genes en `features.tsv` debe coincidir con la primera dimensión en `matrix.mtx`, y el número de códigos de barras en `barcodes.tsv` con la segunda dimensión. Scanpy suele detectar estas inconsistencias.
*   **Recursos del servidor (si está desplegada):** Si la aplicación está desplegada en un servidor con recursos limitados, cargar datasets muy grandes podría exceder la memoria disponible.

**P2: El "Pipeline Principal" falla o se queda bloqueado.**
*   **Datos de entrada de muy baja calidad:** Si los datos tienen muy pocas células o genes después de los filtros iniciales, algunos pasos posteriores (como HVG o PCA) podrían fallar. Revisa tus parámetros de filtrado QC.
*   **Parámetros inadecuados:** Valores extremos en los parámetros del pipeline (ej: pedir 0 HVGs, o un número de PCs mayor al número de células/genes) pueden causar errores.
*   **Memoria insuficiente:** El análisis de scRNA-seq puede consumir mucha memoria, especialmente con datasets grandes. Si se ejecuta localmente, asegúrate de que tu máquina tiene suficientes recursos. Si está desplegada, el servidor podría estar limitando la memoria.
*   **Errores específicos de Scanpy/Anndata:** A veces, pueden surgir errores internos de las bibliotecas. El traceback (registro de error) que muestra la aplicación en rojo puede dar pistas. Si el error no es claro, intenta buscar el mensaje de error específico de Scanpy en internet.

**P3: No se encuentran Genes Altamente Variables (HVGs) o el número es muy bajo.**
*   **Datos muy homogéneos:** Si hay poca variabilidad biológica en tus datos (o es enmascarada por mucho ruido técnico), es posible que se encuentren pocos HVGs.
*   **Parámetros de filtrado QC muy estrictos:** Si has eliminado demasiadas células o genes, podrías estar perdiendo la señal necesaria para detectar HVGs.
*   **Método de detección de HVG:** La aplicación usa `flavor='seurat_v3'`. Para algunos datasets, otros métodos o parámetros dentro de `sc.pp.highly_variable_genes` podrían ser más adecuados (esto requeriría modificar el código).

**P4: El UMAP no muestra una buena separación de clústeres o tiene una forma extraña.**
*   **Número de PCs (`Nº PCs`):** Este es un parámetro crítico. Pocos PCs pueden llevar a una pérdida de información biológica importante. Demasiados PCs pueden introducir ruido. Experimenta con diferentes valores.
*   **Resolución de Leiden (`Resolución Leiden`):** Si la resolución es muy baja, verás pocos clústeres grandes. Si es muy alta, podrías tener muchos clústeres pequeños, a veces fragmentados. Ajusta este parámetro para obtener una granularidad que tenga sentido biológico.
*   **Calidad de los datos y HVGs:** Si el QC no fue óptimo o los HVGs no capturan bien la estructura, el UMAP lo reflejará.
*   **Efectos de batch no corregidos:** Si tienes muestras de diferentes lotes o condiciones que introducen variabilidad técnica fuerte, esto puede dominar el UMAP. Esta aplicación no implementa corrección de batch explícita.

**P5: El Análisis de Expresión Diferencial (DEA) muestra un error de "insuficientes células/grupos".**
*   **Menos de 3 células en un grupo:** El test de Wilcoxon (usado por defecto) requiere un mínimo de células en cada grupo para poder realizar la comparación (Scanpy internamente suele pedir al menos 3). Si después de seleccionar tus condiciones y, opcionalmente, un clúster específico, uno de los grupos resultantes tiene menos de ~3 células, el DEA fallará para esa comparación.
*   **Condiciones mal definidas:** Asegúrate de haber asignado correctamente las muestras a las condiciones y de haber seleccionado dos condiciones diferentes para la comparación.

**P6: Un gen que sé que debería estar, no aparece en los resultados del Explorador de Genes o en los UMAPs de expresión.**
*   **Nombre del gen incorrecto:** Verifica que el nombre/símbolo del gen que has introducido coincida exactamente (sensible a mayúsculas/minúsculas) con los nombres en tu archivo `features.tsv`. A veces `AKT1` es diferente de `Akt1`.
*   **Gen filtrado:** El gen podría haber sido eliminado durante el paso de filtrado de genes (si se expresó en muy pocas células, según el parámetro `Mínimo células/gen`).
*   **El gen no está en el dataset original:** Puede que el gen no haya sido detectado o anotado en tu experimento inicial.

**P7: Los gráficos tardan mucho en generarse o la aplicación va lenta.**
*   **Tamaño del dataset:** Conjuntos de datos con muchas células (>50,000-100,000) pueden hacer que los cálculos y la generación de gráficos sean lentos, especialmente si se ejecuta en un ordenador personal con recursos limitados.
*   **Operaciones de graficación de Matplotlib/Scanpy:** Algunos plots, especialmente los dotplots o violin plots con muchos grupos o características, pueden ser intensivos.

**P8: ¿Cómo puedo guardar todos los resultados de una vez?**
*   La mejor manera de guardar un estado completo del análisis es descargar el **`AnnData Procesado (.h5ad)`** desde la barra lateral. Este archivo contiene los datos normalizados, metadatos, reducciones dimensionales, clusters, etc., y puede ser cargado posteriormente en Python con Scanpy (`sc.read_h5ad("archivo.h5ad")`).
*   Las tablas (marcadores, DEA) se pueden descargar como CSV.
*   Los gráficos individuales se pueden descargar como PNG o HTML.

## 6. Documentación Técnica (para Desarrolladores)

Esta sección proporciona una visión general de la estructura del código, las dependencias clave y el flujo de datos dentro de la aplicación, dirigida a desarrolladores o usuarios avanzados que deseen comprender o modificar el script.

### 6.1. Estructura del Proyecto y Script Principal

El proyecto consiste en un único script de Python (ej: `app_scrna_streamlit.py`) que utiliza la biblioteca Streamlit para generar la interfaz de usuario y las bibliotecas Scanpy y AnnData para realizar los análisis de scRNA-seq.

El script se puede dividir conceptualmente en las siguientes secciones principales:

1.  **Importaciones y Configuración Inicial:**
    *   Importación de las bibliotecas necesarias (`streamlit`, `scanpy`, `anndata`, `matplotlib`, `pandas`, `os`, `tempfile`, `io`, `traceback`, `plotly.express`).
    *   Configuración de la página de Streamlit (`st.set_page_config`).
    *   Título de la aplicación.

2.  **Funciones Auxiliares:**
    *   `load_10x_data()`: Encapsula la lógica para leer datos 10x Genomics desde archivos subidos. Utiliza `tempfile.TemporaryDirectory` para guardar temporalmente los archivos (ya que `sc.read_10x_mtx` espera una ruta de directorio) y luego los lee. Asigna el nombre de la muestra a `adata.obs['sample']`.
    *   `fig_to_bytes()`: Convierte una figura de Matplotlib a un stream de bytes para facilitar su descarga a través de Streamlit.

3.  **Inicialización Centralizada de `st.session_state`:**
    *   Se define un diccionario `default_values` con todas las claves y sus valores por defecto que se usarán en `st.session_state`.
    *   Un bucle asegura que cada clave exista en `st.session_state`, asignándole su valor por defecto si no está ya presente. Esto es crucial para mantener el estado de la aplicación entre interacciones y evitar errores de `KeyError`.

4.  **Barra Lateral (Sidebar - `st.sidebar`):**
    *   Contiene la mayoría de los controles de entrada del usuario.
    *   **Sección 1: Carga de Datos:**
        *   Widget `st.number_input` para el número de muestras.
        *   Bucle para generar dinámicamente los `st.text_input` (nombre de muestra) y `st.file_uploader` (para matrix, features, barcodes) para cada muestra. Los archivos y nombres se almacenan en `st.session_state.sample_files`.
    *   **Sección 2: Parámetros de Pipeline Principal:**
        *   Widgets `st.slider`, `st.text_input` para configurar los parámetros de QC, HVG, PCA, clustering. Los valores se almacenan directamente en `st.session_state` (ej: `st.session_state.min_genes_val`).
    *   **Sección 3: Análisis de Expresión Diferencial (DEA):**
        *   Esta sección se renderiza condicionalmente (`if st.session_state.analysis_done and st.session_state.adata_processed is not None`).
        *   Widgets para asignar condiciones, seleccionar grupos de comparación, y definir parámetros del DEA.
        *   Botón `st.button` para ejecutar el DEA, que llama a la lógica de `sc.tl.rank_genes_groups`.
    *   **Botones de Acción Principales:**
        *   `Cargar y Concatenar Datos`: Ejecuta la carga usando `load_10x_data` y `ad.concat`. Actualiza `st.session_state.adata_combined_raw`.
        *   `Ejecutar Pipeline Principal`: Ejecuta la secuencia completa de análisis de Scanpy. Actualiza `st.session_state.adata_processed`, `st.session_state.marker_genes_df`, etc.
    *   Botón de descarga para el `AnnData` procesado.

5.  **Sección de Resultados (Panel Principal):**
    *   Encabezado y el `st.text_area` para el "Explorador de Expresión Génica", cuyo valor se guarda en `st.session_state.gene_explorer_input`.
    *   Renderizado condicional principal (`if st.session_state.analysis_done and st.session_state.adata_processed is not None`).
    *   Uso de `st.tabs` para organizar los diferentes tipos de resultados.
    *   **Dentro de cada pestaña:**
        *   Lógica para generar los plots usando funciones de `sc.pl` (ej: `sc.pl.umap`, `sc.pl.dotplot`, `sc.pl.violin`). Las figuras de Matplotlib se muestran con `st.pyplot()`.
        *   Las tablas de Pandas (marcadores, DEA) se muestran con `st.dataframe()`.
        *   El Volcano Plot se genera con `plotly.express` y se muestra con `st.plotly_chart()`.
        *   Botones `st.download_button` para descargar figuras y tablas.
    *   Manejo de errores con `try-except` alrededor de las operaciones de ploteo y análisis para evitar que la aplicación se bloquee por completo.

### 6.2. Dependencias Principales

El archivo `requirements.txt` (si se proporciona) listaría las dependencias exactas. Las bibliotecas fundamentales son:

*   `streamlit`: Para la interfaz web.
*   `scanpy`: Para el análisis de scRNA-seq.
*   `anndata`: Estructura de datos subyacente a Scanpy.
*   `matplotlib`: Para la generación de gráficos estáticos (usado por Scanpy).
*   `pandas`: Para la manipulación de datos tabulares.
*   `plotly` (específicamente `plotly.express`): Para gráficos interactivos como el Volcano Plot.
*   `numpy` (implícita, dependencia de Scanpy/Pandas).
*   `scipy` (implícita, dependencia de Scanpy).

Se recomienda instalar estas dependencias en un entorno virtual de Python.
Ej: `pip install streamlit scanpy anndata matplotlib pandas plotly`

### 6.3. Flujo de Datos (Objetos AnnData Clave)

La aplicación gestiona varios estados del objeto AnnData a través de `st.session_state`:

1.  **`st.session_state.sample_files` (dict):** Almacena los objetos de archivo subidos y los nombres de muestra antes de la concatenación.
2.  **`st.session_state.adata_combined_raw` (AnnData):**
    *   Creado al pulsar "Cargar y Concatenar Datos".
    *   Contiene los datos crudos concatenados de todas las muestras.
    *   Las células tienen `adata.obs['sample']` para indicar su origen.
    *   Sirve como punto de partida para el pipeline principal.
3.  **`adata_pipeline` (AnnData, variable local dentro de la función del pipeline):**
    *   Una copia de `adata_combined_raw` sobre la cual se ejecutan los pasos de QC, normalización y cálculo de HVGs (sobre todos los genes que pasan QC).
4.  **`st.session_state.adata_hvg_filtered_intermediate` (AnnData):**
    *   Creado durante el "Pipeline Principal".
    *   Una subselección de `adata_pipeline` que contiene *solo* los genes altamente variables (HVGs) en sus columnas (`.var`) y las células que pasaron el QC en sus filas (`.obs`).
    *   Los datos en `.X` de este objeto son los que se escalan y sobre los que se ejecuta PCA.
    *   Los resultados de PCA, cálculo de vecinos y UMAP se calculan sobre este objeto.
5.  **`st.session_state.adata_processed` (AnnData):**
    *   El objeto AnnData final que se utiliza para la mayoría de las visualizaciones.
    *   Se basa en `adata_pipeline` (es decir, contiene todos los genes que pasaron el QC inicial, no solo los HVGs).
    *   Los resultados clave del análisis sobre `adata_hvg_filtered_intermediate` se transfieren a este objeto:
        *   `adata_processed.obsm['X_umap']` se copia de `adata_hvg_filtered_intermediate.obsm['X_umap']`.
        *   `adata_processed.obs['leiden_clusters']` se copia de `adata_hvg_filtered_intermediate.obs['leiden_clusters']`.
    *   Los genes marcadores de clúster (`sc.tl.rank_genes_groups`) se calculan sobre este `adata_processed` usando los datos normalizados y logaritmizados (pero no necesariamente escalados o filtrados por HVG para este paso específico, dependiendo de `use_raw`).
    *   Si se realiza DEA, se hace una copia de este objeto para añadir la columna `condition` y realizar el `rank_genes_groups`.

### 6.4. Puntos de Extensión o Modificación Sugeridos

*   **Corrección de Efectos de Batch:** Integrar métodos de Scanpy como `sc.external.pp.bbknn` o `sc.external.pp.harmony_integrate` después de la carga de datos y antes del cálculo de vecinos/UMAP. Esto requeriría ajustar el flujo de datos y la interfaz.
*   **Más Opciones de Normalización/Escalado/HVG:** Exponer más parámetros o métodos alternativos de Scanpy en la interfaz.
*   **Visualizaciones Avanzadas:** Añadir plots como "elbow plot" para la selección de PCs, trayectorias de pseudotiempo (si es aplicable), o heatmaps más personalizables.
*   **Anotación Celular:** Integrar herramientas o estrategias para la anotación de tipos celulares (ej: carga de listas de marcadores, métodos automatizados si existen wrappers en Scanpy).
*   **Interactividad Mejorada en Plots:** Donde sea posible, usar Plotly para más gráficos (ej: UMAPs interactivos) en lugar de Matplotlib estático.
*   **Guardar/Cargar Estado de la Sesión:** Permitir al usuario guardar la configuración de parámetros y los resultados intermedios para reanudar un análisis más tarde (Streamlit no ofrece esto de forma nativa fácilmente para objetos complejos como AnnData).

## 7. Glosario de Términos

Aquí se definen algunos términos clave utilizados en la aplicación y en el campo del análisis de single-cell RNA-seq.

*   **scRNA-seq (Single-Cell RNA Sequencing / Secuenciación de ARN de Célula Única):**
    Técnica que permite medir los niveles de expresión génica (el transcriptoma) de miles de células individuales de forma simultánea.

*   **AnnData (Annotated Data / Datos Anotados):**
    Es la estructura de datos central utilizada por Scanpy (y otras herramientas) para almacenar datos de scRNA-seq. Un objeto AnnData contiene la matriz de expresión, anotaciones para células (observaciones, `.obs`), anotaciones para genes (variables, `.var`), y resultados de análisis dimensional (en `.obsm`), entre otros.

*   **10x Genomics:**
    Una compañía que proporciona instrumentación y reactivos populares para experimentos de scRNA-seq. El formato de datos de salida de su plataforma Cell Ranger (archivos `matrix.mtx`, `features.tsv`, `barcodes.tsv`) es un estándar de facto.

*   **QC (Quality Control / Control de Calidad):**
    Proceso de identificar y eliminar células o genes de baja calidad de los datos. Las métricas comunes incluyen el número de genes detectados por célula, el número total de cuentas (UMIs) por célula, y el porcentaje de cuentas de genes mitocondriales.

*   **UMI (Unique Molecular Identifier / Identificador Molecular Único):**
    Una etiqueta de secuencia corta que se añade a las moléculas de ARN antes de la amplificación por PCR. Permite contar el número real de moléculas de transcritos originales, corrigiendo los sesgos de amplificación. Las "cuentas" en scRNA-seq suelen referirse a cuentas de UMI.

*   **Normalización:**
    Proceso para ajustar las cuentas de expresión crudas para eliminar diferencias técnicas entre células (ej: diferencias en la eficiencia de captura de ARN o profundidad de secuenciación), permitiendo comparaciones más justas de los niveles de expresión. Un método común es la normalización por tamaño de librería (ej: escalar cada célula a 10,000 cuentas totales).

*   **Log Transformación (Transformación Logarítmica, ej: `log1p`):**
    Aplicación de una función logarítmica (a menudo `log(x+1)`) a los datos normalizados. Ayuda a estabilizar la varianza y a hacer que los datos se asemejen más a una distribución normal, lo cual es beneficioso para algunos algoritmos downstream.

*   **HVG (Highly Variable Genes / Genes Altamente Variables):**
    Genes que muestran una variabilidad en su expresión entre células mayor de la esperada por azar. Se presume que esta variabilidad refleja diferencias biológicas y son los más informativos para distinguir tipos celulares o estados.

*   **PCA (Principal Component Analysis / Análisis de Componentes Principales):**
    Técnica de reducción de dimensionalidad que transforma los datos (generalmente los HVGs) en un nuevo conjunto de variables no correlacionadas llamadas componentes principales (PCs), ordenadas por la cantidad de varianza que explican. Se utiliza para reducir el ruido y la dimensionalidad antes de pasos como UMAP o clustering.

*   **UMAP (Uniform Manifold Approximation and Projection):**
    Algoritmo de reducción de dimensionalidad no lineal que se utiliza comúnmente para visualizar datos de scRNA-seq en 2D o 3D. Intenta preservar tanto la estructura global como la local de los datos.

*   **Clustering (Agrupamiento):**
    Proceso de agrupar células similares basándose en sus perfiles de expresión génica. El objetivo es identificar poblaciones celulares distintas.
    *   **Leiden Clustering:** Un algoritmo popular para la detección de comunidades en grafos, comúnmente aplicado al grafo de vecindad de células en scRNA-seq. La `Resolución Leiden` controla la granularidad de los clústeres.

*   **Genes Marcadores (Marker Genes):**
    Genes que se expresan de forma diferencial y preferente en un clúster celular (o tipo celular) en comparación con otros. Son cruciales para identificar y caracterizar los clústeres.

*   **DEA (Differential Expression Analysis / Análisis de Expresión Diferencial):**
    Análisis estadístico para identificar genes cuya expresión cambia significativamente entre dos o más condiciones o grupos de células (ej: control vs. tratado, o tipo celular A vs. tipo celular B).

*   **Log2FC (Log2 Fold Change / Cambio Logarítmico en Base 2):**
    Una medida de la magnitud del cambio en la expresión génica. Un Log2FC de 1 significa un aumento de 2 veces en la expresión; un Log2FC de -1 significa una disminución de 2 veces (la mitad de expresión).

*   **P-valor:**
    En DEA, la probabilidad de observar un cambio de expresión tan grande o mayor que el medido, si no hubiera una diferencia real entre los grupos (bajo la hipótesis nula). Valores pequeños sugieren que el cambio observado es improbable por azar.

*   **P-valor Ajustado (Adjusted P-value / P-adj, q-value):**
    P-valor corregido para tener en cuenta las múltiples pruebas estadísticas realizadas (una por cada gen). Ayuda a controlar la tasa de falsos descubrimientos. Es el valor que se suele usar para determinar la significancia estadística en DEA.

*   **Dot Plot (Diagrama de Puntos):**
    Un tipo de gráfico utilizado en scRNA-seq para visualizar la expresión de múltiples genes a través de múltiples grupos (clústeres o condiciones). Típicamente, el color del punto representa el nivel promedio de expresión y el tamaño del punto representa el porcentaje de células en el grupo que expresan el gen.

*   **Violin Plot (Diagrama de Violín):**
    Un gráfico que combina características de un box plot y un gráfico de densidad de kernel. Muestra la distribución de datos numéricos, útil para comparar distribuciones de expresión génica entre grupos.


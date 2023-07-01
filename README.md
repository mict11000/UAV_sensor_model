# Tutorial para modelar un sensor pasivo (cámaras) y activo (Lidar) en un RPAS 

Este repositorio contiene el código en Python de un módulo que implementa clases para modelar un sensor pasivo es decir una cámara digital, permitiendo obtener algunas de sus propiedades, cálculo de las coordenadas de los vértices de una imagen aérea y su proyección según su posición y orientación, la huella en el terreno (footprint) como un vector que se puede visualizar en un SIG, por ultimo permite aplicar la georreferenciación directa, se utilizó las librerías fundamentales de Python como GDAL Y Numpy.

# Instalación y creación del ambiente

Descargar e instalar la versión de Anaconda para tu sistema operativo desde su página web.

Ejecutar en Anaconda Prompt o el terminal, las siguientes líneas de código:

Crear un nuevo ambiente:

```
conda create --name osgeo_env -c conda-forge python=3.10
```

Instalar las librerías necesarias:
```
conda install -c conda-forge mamba
```

```
mamba install -c conda-forge gdal
```

```
mamba install -c conda-forge geopandas
```

```
mamba install -c conda-forge jupyterlab
```

Ejecutar jupyter lab:

 ```
jupyter lab
```






## Referencias

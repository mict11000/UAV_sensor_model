# Tutorial para modelar un sensor pasivo (cámaras) y activo (Lidar) en un RPAS 

Este repositorio contiene el código en Python de un módulo que implementa diferente clases; una para modelar un sensor pasivo es decir una cámara digital, permitiendo obtener algunas de sus propiedades, cálculo de las coordenadas de los vértices de una imagen aérea y su proyección según su posición y orientación, la huella en el terreno (footprint) como un vector tipo polígono que se puede visualizar en un SIG, por ultimo permite aplicar la georreferenciación directa. Adicionalmente existe otras dos clases aun en desarrollo, una para modelar la algunas propiedades del sistema Lidar, y para obtener su huella, la ruta de vuelo como un vector tipo línea a partir de una lista de coordenadas y el cambio de dirección entre lineas de vuelo mediante una semicircunferencia generada a partir de un polígono regular inscrito en una circunferencia, se utilizó las librerías fundamentales de Python como GDAL Y Numpy.

Para poder demostrar la funcionalidad del módulo se provee un notebook que puede ser ejecutado en el siguiente ambiente: 

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
jupyter lab --browser=firefox
```


## Referencias

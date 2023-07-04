# Tutorial para modelar un sensor pasivo (cámaras) y activo (Lidar) en un RPAS 

Este repositorio contiene el código en Python de un módulo que implementa diferente clases: una para modelar un sensor pasivo es decir una cámara digital, permitiendo obtener algunas de sus propiedades, las coordenadas de los vértices de una imagen aérea y su proyección según su posición y orientación, la huella en el terreno (footprint) como un vector tipo polígono que se puede visualizar en un SIG, por ultimo permite aplicar la georreferenciación directa a la imagen, esto algoritmos estan basados en supuestos que simplifican la demostración de las ecuaciones. Adicionalmente existen otras dos clases aun en desarrollo, una para modelar algunas propiedades del sistema Lidar, y obtener su huella. La ultima clase ayuda con el diseño de la ruta de vuelo como un vector tipo línea a partir de una lista de coordenadas y el cambio de dirección entre lineas de vuelo mediante una semicircunferencia generada a partir de un polígono regular inscrito en una circunferencia, se utilizó las librerías fundamentales de Python como GDAL, Numpy, GeoPandas, Folium y Jupyter Lab.

Para poder demostrar la funcionalidad del módulo se provee un notebook que puede ser ejecutado en el siguiente ambiente: 

# Instalación y creación del ambiente

Descargar e instalar la versión de Anaconda para tu sistema operativo desde su página web.

Ejecutar en Anaconda Prompt o el terminal, las siguientes líneas de código:

Crear un nuevo ambiente, con el administrador de paquetes mamba:

```
conda create -n osgeo_env -c conda-forge python=3.10 mamba
```

Activar el ambiente previamente creado:

```
conda activate osgeo_env
```

Instalar las librerías necesarias:

```
mamba install -c conda-forge gdal geopandas folium jupyterlab
```

Ejecutar jupyter lab, con el navegaor de nuestra preferencia:

 ```
jupyter lab --browser=firefox
```


## Referencias

# Modelamiento del Balance de Energía - Glaciar Juncal Norte

Este repositorio contiene el código fuente desarrollado en Google Earth Engine (GEE) para el modelamiento de la radiación de onda corta (basado en Zhang et al., 2015) y el cálculo de flujos turbulentos (Schaefer, 2020). 

## Contenido
* `modelo_balance_juncal.js`: Script principal de GEE que realiza la corrección topográfica, cálculo de albedo satelital y balance energético (Qm).

## Uso
El código está diseñado para ejecutarse en el Code Editor de Google Earth Engine. Requiere acceso a la colección Sentinel-2 (Harmonized) y un DEM de 10m de resolución.

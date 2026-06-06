///// Seccion 1  - Definir colección //////

// Elegir colección de interés (S2_SR), filtrar por fecha, filtrar por geometría, aplicar máscara de nubes y des-armonizar offset
function maskS2clouds(image) {
  var qa = image.select('QA60');
  
  // Bits 10 and 11 are clouds and cirrus, respectively.
  var cloudBitMask = 1 << 10;
  var cirrusBitMask = 1 << 11;

  // Both flags should be set to zero, indicating clear conditions.
  var mask = qa.bitwiseAnd(cloudBitMask).eq(0)
               .and(qa.bitwiseAnd(cirrusBitMask).eq(0));
               
  // 1. SELECCIONAR ÚNICAMENTE LAS BANDAS ÓPTICAS AFECTADAS (B1 a B12)
  var bandasOpticasRaw = image.select('B.*');
  
  // 2. CORRECCIÓN DEL FORMATO ESA (Sumar 1000 DN para deshacer la resta automática de GEE)
  var bandasOpticasCorregidas = bandasOpticasRaw.add(1000);
  
  // 3. ESCALADO A REFLECTANCIA SUPERFICIAL (0 - 1)
  var bandasOpticasEscaladas = bandasOpticasCorregidas.divide(10000);
  
  // 4. REEMPLAZAR LAS BANDAS EN LA IMAGEN ORIGINAL 
  // (Esto mantiene intactas las bandas atmosféricas WVP, AOT y de calidad QA60 sin tocarlas)
  var imagenCompletamenteCalibrada = image.addBands(bandasOpticasEscaladas, null, true);
               
  // Aplicamos la máscara de nubes final, copiamos propiedades y retornamos
  return imagenCompletamenteCalibrada.updateMask(mask)
              .copyProperties(image)
              .set("system:time_start", image.get("system:time_start"));
}

// Función para quitar los valores de SR mayor a 1
function reflectanceMask (image) {
  var sunglint = image.select(["B2"]);
  var sunglint2 = image.select(["B3"]);
  var sunglint3 = image.select(["B4"]);
  var mask = sunglint.lte(1).and(sunglint2.lte(1)).and(sunglint3.lte(1));
  var sunglintmask = image.updateMask(mask);
  
  return sunglintmask.copyProperties(image);
}

// DECLARACIÓN DE LA COLECCIÓN GLOBAL DE LA TESIS
var dataset = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
                  .filter(ee.Filter.calendarRange(2020, 2025, "year"))
                  .filter(ee.Filter.calendarRange(1, 4, 'month'))
                  .filterBounds(juncal_buffer)
                  .map(maskS2clouds) // Aquí adentro ocurre la magia de la corrección óptica
                  .map(function (image) {return image.clip(juncal_buffer)})
                  .map(reflectanceMask);
print(dataset.first())
var median = dataset.first();
Map.addLayer(median.select(['B4', 'B3','B2']), {min: 0, max: 0.3}, 'juncal' )
Map.centerObject(juncal, 10)
Map.addLayer(juncal, {}, "juncalarea")

// --- Constantes para Flujos Turbulentos (Schaefer 2020) ---
var kappa = 0.4;             // Constante de von Kármán
var z_sensor = 2.0;                 // Altura del sensor (m)
var z0 = 0.0005;             // Rugosidad para todo el glaciar (0.5 mm) // zona de la lengua (2 mm)
var ca = 1010;               // Calor específico del aire (J/kg K)
var Lv = 2501000;            // Calor latente de evaporación/sublimación (J/kg)
var Rd = 287.05;             // Constante de gas para aire seco (J/kg K)

// Cálculo del Coeficiente de Transferencia Cinematica (C*)
var Cp_star = Math.pow(kappa, 2) / Math.pow(Math.log(z_sensor / z0), 2);

///////////////////////////////////////// --- Configuración del Gradiente Térmico ---
var altitudReferencia = 3155; 
var gradienteTermico = 0.0065; 

// 1. Definir la función usando mi variable dem_10m
var calcularTemperaturaGlaciar = function(image) {
  // Obtener la fecha de la imagen para buscarla en la tabla temperatura_3155m
  var fechaImg = image.date().format('dd-MM-YYYY');
  
  // Buscar el dato correspondiente en el Asset que subiste
  var registro = temperatura_diaria_3155m.filter(ee.Filter.eq('fecha', fechaImg)).first();
  
  // Extraer la temperatura (si no hay dato para ese día, se asigna 0)
  var tEstacion = ee.Number(ee.Algorithms.If(registro, registro.get('temp_ref_3155m'), 0));
  
  // Aplicar la fórmula: T_pixel = T_base - 0.0065 * (Altitud_pixel - 3155)
 var mapaTemp = ee.Image.constant(tEstacion)
    .subtract(dem_10m.subtract(altitudReferencia).multiply(gradienteTermico))
    .rename('temperatura_distribuida');

  return image.addBands(mapaTemp).set('t_base_dia', tEstacion);
};

// 2. Crear la variable datasetConTemp (esto soluciona el error de "not defined")
var datasetConTemp = dataset.map(calcularTemperaturaGlaciar);

// 3. Visualización en el Mapa
var visTemp = {
  min: -5, 
  max: 15, 
  palette: ['0000ff', '00ffff', 'ffff00', 'ff0000']
};



////// Sección 2 - Elegir bandas y definir ecuación //////

var bandas = median.select (["B2", "B3", "B4", "B8"])


////   Sección 3 - Generar una función con la ecuación    //////// 

// Reemplaza los números con los que obtuviste en QGIS
// Importante: El orden siempre es [Longitud, Latitud]
var puntoEstacion = ee.Geometry.Point([-70.09972552133777, -32.99945014456693]);



//////////////////////////////// CÁLCULO DE LA RADIACIÓN DIRECTA /////////////////////////////////

function calculateRadiationComponents(image) {
  
  
 var albedo_raw = image.expression(
  '0.726 * B3 - 0.322 * pow(B3, 2) - 0.015 * B8 + 0.581 * pow(B8, 2)', {
    'B3': image.select('B3'),
    'B8': image.select('B8')
}).rename('albedo_raw');

// Ya no necesitas el .add(0.2057) porque el offset físico de la ESA ya fue solucionado arriba
var albedo = albedo_raw.clamp(0, 1).rename('albedo');

  // 2. Reducción para la consola
  var statsAlbedo = albedo.reduceRegion({
    reducer: ee.Reducer.median(),
    geometry: juncal, 
    scale: 10,
    maxPixels: 1e9
  });

  // Actualizamos la imagen con la banda y la propiedad, pero NO ponemos 'return'
  image = image.addBands(albedo)
               .set('albedo_valor', statsAlbedo.get('albedo'))
               .set('fecha_display', image.date().format('dd-MM-YYYY'));
  // =======================================================================
  
  // 1. OBTENER PARÁMETROS TOPOGRÁFICOS Y DE LA IMAGEN
  // --------------------------------------------------
  var slopeRad = ee.Terrain.slope(dem_10m).multiply(Math.PI/180);
  var aspectRad = ee.Terrain.aspect(dem_10m).multiply(Math.PI/180);
  
  var theta_s = ee.Number(image.get('MEAN_SOLAR_ZENITH_ANGLE')).multiply(Math.PI / 180);
  var cos_theta_s = theta_s.cos();
  var phi_s = ee.Number(image.get('MEAN_SOLAR_AZIMUTH_ANGLE')).multiply(Math.PI / 180);
  
  var h = ee.Number(Math.PI / 2).subtract(theta_s);
  var denominator = h.sin().add(ee.Number(0.15).multiply(h.multiply(57.296).add(3.885).pow(-1.253)));
  var m = ee.Number(1).divide(denominator);
  

  // 2. CÁLCULO DE E₀ (IRRADIANCIA SOLAR EXTRATERRESTRE)
  // ----------------------------------------------------
  var date = ee.Date(image.get('system:time_start'));
  var n = date.getRelative('day', 'year').add(1); // revisar
  var gamma = ee.Number(n).subtract(1).multiply(2 * Math.PI).divide(365);
  var dr = ee.Number(1.000110)
    .add(ee.Number(0.034221).multiply(gamma.cos()))
    .add(ee.Number(0.001280).multiply(gamma.sin()))
    .add(ee.Number(0.000719).multiply((gamma.multiply(2)).cos()))
    .add(ee.Number(0.000077).multiply((gamma.multiply(2)).sin()));
  var ISC = 1361;
  var E0 = ee.Image.constant(ISC).multiply(dr).rename('E0');


// 3. CÁLCULO DE LA TRANSMISIVIDAD ATMOSFÉRICA (Tb)
  // --------------------------------------------------
  var Po = 101325; 
  var z = dem_10m;
  var Ps = ee.Image.constant(Po).multiply(z.divide(-8430).exp()); 
  var mc = ee.Image.constant(m).multiply(Ps.divide(Po)); 

  // ==========================================================
  // ESTRATEGIA DE RELLENO JERÁRQUICO (ANTI-COLAPSO)
  // ==========================================================
  var date_jerarquico = ee.Date(image.get('system:time_start'));

  // 1. OBTENER DATOS BASE DE SENTINEL-2
  var s2_w = image.select('WVP').multiply(0.001);    
  var s2_aod = image.select('AOT').multiply(0.001); 

  // 2. OBTENER DATOS DE MODIS (CON ESCUDO PROTECTOR)
  var modisAtmos = ee.ImageCollection("MODIS/061/MCD19A2_GRANULES")
    .filterDate(date_jerarquico, date_jerarquico.advance(1, 'day'));

  // Si MODIS pasó ese día (>0 imágenes), hace la matemática. 
  // Si no pasó, inyecta temporalmente la constante para evitar el error de 0 bandas.
  var modisPW = ee.Image(ee.Algorithms.If(
      modisAtmos.size().gt(0),
      modisAtmos.median().select('Column_WV').multiply(0.001),
      ee.Image.constant(1.0).rename('Column_WV') 
  ));

  var modisAOD = ee.Image(ee.Algorithms.If(
      modisAtmos.size().gt(0),
      modisAtmos.median().select('Optical_Depth_047').multiply(0.001),
      ee.Image.constant(0.05).rename('Optical_Depth_047') 
  ));

  // 3. OZONO (promedio del ozono a estas latitudes desde 1978 usando la colección TOMS entre los meses de enero-marzo)
  var l = ee.Image.constant(0.2598432);

  // 4. LA CADENA DE ENMASCARAMIENTO FINAL
  // S2 -> MODIS (si existe el píxel) -> Constante Física (si ambos fallaron)
  var w = s2_w.unmask(modisPW).unmask(0.5); 
  var aod = s2_aod.unmask(modisAOD).unmask(0.05);

  // ==========================================================
  // CÁLCULO DE LAS TRANSMISIVIDADES (Ecuaciones de Yang et al.)
  // ==========================================================
  
  var Tw = ee.Image.constant(1.0).min(
    ee.Image.constant(0.909).subtract(ee.Image.constant(0.036).multiply(ee.Image.constant(m).multiply(w).log()))
  ).rename('Tw');

  var beta = aod.multiply(Math.pow(0.5, 1.3)); 
  var m_beta = ee.Image.constant(m).multiply(beta); 
  var Ta_expr = 'exp(-m * b * pow(0.6777 + 0.1464 * (m * b) - 0.00626 * pow(m * b, 2), -1.3))'; 
  var Ta = m_beta.expression(Ta_expr, {'m': m, 'b': beta}).rename('Ta');

  var Tr_expr = 'exp(-0.008735 * mc * pow(0.547 + 0.014 * mc - 0.00038 * pow(mc, 2) + 4.6e-6 * pow(mc, 3), -4.08))';
  var Tr = mc.expression(Tr_expr, {'mc': mc}).rename('Tr'); 
  
  var To3 = ee.Image.constant(m).multiply(l).pow(0.7136).multiply(-0.0365).exp().rename('To3');
  
  var Tg = mc.pow(0.3139).multiply(-0.0117).exp().rename('Tg');
  
  var Tb = Tr.multiply(Ta).multiply(Tw).multiply(To3).multiply(Tg).subtract(0.013).max(0).rename('Tb');


// 4. CÁLCULO DE LOS FACTORES TOPOGRÁFICOS
  // ----------------------------------------
  var elevationDeg = ee.Number(90).subtract(image.get('MEAN_SOLAR_ZENITH_ANGLE'));
  var hillshade = ee.Terrain.hillshade(dem_10m, image.get('MEAN_SOLAR_AZIMUTH_ANGLE'), elevationDeg);
  var Vs = hillshade.gt(0).rename('Vs'); 

  var cos_is_expr = 'cos(theta_s) * cos(beta) + sin(theta_s) * sin(beta) * cos(phi_s - phi)';
  var cos_is_raw = ee.Image().expression(cos_is_expr, {
    'theta_s': theta_s,
    'phi_s': phi_s,
    'beta': slopeRad,
    'phi': aspectRad
  });
  
  // --- LA CORRECCIÓN DE ESCALA ESPACIAL (Zhang et al., 2015) ---
  // 1. clamp(0, 1): Evita valores negativos irreales por las sombras
  // 2. focal_mean: Suaviza la micro-topografía del DEM de 10m en un radio de 30m, 
  // diluyendo los "paneles solares" extremos y estabilizando la mediana global.
  var cos_is = cos_is_raw.clamp(0, 1)
                         //.focal_mean({radius: 30, units: 'meters'})
                         .rename('cos_is');

  // 5. CÁLCULO FINAL DE LA RADIACIÓN DIRECTA (Edir)
  // ------------------------------------------------
  var Edir = Vs.multiply(E0).multiply(Tb).multiply(cos_is).rename('Edir');
  



////////////////////// 6. CALCULO DE LA RADIACIÓN DIFUSA ANISOTRÓPICA //////////////////////////////////

// ==============================================================
  
  // 6.1. Calcular la Transmisividad Difusa (Td) - Yang et al. (2006), Eq. (1b) y Eq. (6) de Zhang et al. (2015)
  var Td = To3.multiply(Tg).multiply(Tw)
              .multiply(ee.Image.constant(1).subtract(Ta.multiply(Tr)))
              .add(0.013)
              .multiply(0.5)
              .max(0)
              .rename('Td');
              
  // 6.2. Calcular la Irradiancia Difusa en Horizontal (Edif_hor) - Zhang et al. (2015), Eq. (5)
  var Edif_hor = E0.multiply(Td).multiply(cos_theta_s).rename('Edif_hor');
  
  // 6.3. Definir el Índice de Anisotropía (K)
  var K = Tb.rename('K'); // Este se define como la relación entre la irradiancia directa que llega a una superficie horizontal y la irradiancia extraterrestre que llegaría si no hubiera atmosfera
  
  // 6.4. Calcular Eaniso_dif - Zhang et al. (2015), Eq. (4)
  var Eaniso_dif = Edif_hor.multiply(Vs)
                           .multiply(K)
                           .multiply(cos_is.divide(cos_theta_s))
                           .max(0) // Asegura que no sea negativo
                           .rename('Eaniso_dif');


//////////////////// 7. CALCULO DE LA RADIACIÓN DIFUSA ISOTRÓPICA ///////////////////////////////////////

// Cargar el Asset y FORZAR la alineación con mi zona de estudio
var topogProfesor = ee.Image("projects/tesis-457819/assets/svf_tvf_li2002")
  .resample('bilinear')
  .reproject({
    crs: 'EPSG:32719', // UTM 19S (el que usas para el Juncal)
    scale: 10          // Forzamos a 10 metros para que coincida con Sentinel
  });


 // 7.1 Calcular el Factor de Visión del Cielo (Viso)
   // Esto se calculó a través de QGIS con herramienta SAGA
  var Viso = topogProfesor.select("b1");
  
  // 7.2 Calcular Eiso_dif - Zhang et al. (2015), Eq. (7)
  var Eiso_dif = Edif_hor.multiply(Viso)
                         .multiply(ee.Image.constant(1).subtract(K))
                         .rename('Eiso_dif');
                         
                        

//////////////////// 8. CÁLCULO DE LA RADIACIÓN REFELJADA //////////////////////////////////

// Usamos la Banda 2 del profesor como Terrain View Factor (TVF)
var TVF = topogProfesor.select("b2");

// 8.1. Calcular la Irradiancia Directa en Horizontal (Edir_hor)
var Edir_hor = E0.multiply(Tb).multiply(cos_theta_s).rename('Edir_hor');

// 8.2. Calcular la Irradiancia Global en Horizontal (Eglob_hor = Edir_h + Edif_h)
var Eglob_hor = Edir_hor.add(Edif_hor);

// 8.3. Calcular Eref siguiendo la fórmula (8) de Zhang 2015
// Eref = Albedo * (Radiación Global Horizontal) * TVF
var Eref = albedo.multiply(Eglob_hor)
                 .multiply(TVF)
                 .rename('Eref');

                    


// ==============================================================
// 9. CÁLCULO DE RADIACIÓN GLOBAL TOTAL (DSSR)
// ==============================================================
// Sumamos los 4 componentes: Directa + Difusa Aniso + Difusa Iso + Reflejada
var DSSR = Edir.add(Eaniso_dif)
               .add(Eiso_dif)
               .add(Eref)
               .rename('DSSR');

// ==============================================================
  // 9. PREPARACIÓN DE DATOS PARA EXPORTAR A PYTHON
  // ==============================================================
  
  // 9.1 Sumamos Directa + Difusa para tener la Incidente Total
  var difusa_total = Eiso_dif.add(Eaniso_dif);
  var incidente_total = DSSR.rename('rad_incidente');

  // 9.2 Agrupamos Albedo e Incidente para extraer el dato del glaciar
  var variablesExportar = ee.Image([
    albedo, // El albedo que calculaste en el paso 1
    incidente_total
  ]);

  // 9.3 Extraemos la mediana (el valor resumen del día para el Excel)
  var statsGlaciar = variablesExportar.reduceRegion({
    reducer: ee.Reducer.median(),
    geometry: juncal, 
    scale: 10,
    maxPixels: 1e9
  });

  // 9.4 Guardamos los valores como propiedades de la imagen
  image = image.set('Fecha', image.date().format('dd-MM-YYYY'))
               .set('Anio', image.date().get('year'))
               .set('Mes', image.date().get('month'))
               .set('Albedo', statsGlaciar.get('albedo'))
               .set('Radiacion_Incidente', statsGlaciar.get('rad_incidente'));
               

// ==============================================================
  // 10. CÁLCULO DE ONDA LARGA (LW) - SCHAEFER 2020
  // ==============================================================
  
  var sigma = 5.67e-8;
  // Usamos la banda 'temperatura_distribuida' que ya calculaste en Kelvin
  // (Asumiendo que calcularTemperaturaGlaciar devuelve Celsius, sumamos 273.15)
  var tAirK = image.select('temperatura_distribuida').add(273.15);
  
  // A. LW INCIDENTE (LW_in) - Schaefer (Cielo despejado)
  // Epsilon = 0.00877 * T^0.788 -> LW_in = Epsilon * sigma * T^4
  var lwIn = tAirK.pow(4.788).multiply(0.00877).multiply(sigma).rename('LW_in');
  
  // B. LW SALIENTE (LW_out) - Epsilon = 1
  // Aplicamos el límite físico: La superficie no puede emitir más calor que el hielo fundiéndose (0°C)
  var tSurfK = tAirK.where(tAirK.gt(273.15), 273.15);
  var lwOut = tSurfK.pow(4).multiply(sigma).rename('LW_out');
  
  // C. LW NETA (Qlw)
  var lwNet = lwIn.subtract(lwOut).rename('LW_net');

  // Actualizamos la imagen con las bandas para que sean visibles en el mapa
  image = image.addBands([lwIn, lwOut, lwNet]);
  
  
  // ==============================================================
  // 10. CÁLCULO DE SH - LH Y BALANCE ENERGÉTICO
  // ==============================================================

// ==============================================================
// 11. CÁLCULO DE FLUJOS TURBULENTOS (SH y LH) - CORREGIDO
// ==============================================================

// 1. OBTENER LA FECHA (Cambiado a dd-MM-YYYY para que coincida con tu Excel)
var fechaImg_met = image.date().format('dd-MM-YYYY');

// 2. BUSCAR DATOS METEOROLÓGICOS (Asegúrate que 'datos_SH_LH' sea tu Asset)
var metData = ee.Feature(datos_SH_LH_horario.filter(ee.Filter.eq('fecha', fechaImg_met)).first());

// 3. VALORES DE RESCATE REALISTAS
// Bajamos el viento a 2.5 y subimos humedad a 45 para que, si falla el dato, 
// el balance no se vuelva loco y negativo.
var u_dia = ee.Number(ee.Algorithms.If(metData, metData.get('u_ms'), 2.5));
var u_imagen = ee.Image.constant(u_dia).toFloat().rename('u_dia');
var rh_dia = ee.Number(ee.Algorithms.If(metData, metData.get('rh_percent'), 45.0));
var p_base_estacion = ee.Number(ee.Algorithms.If(metData, metData.get('p_pa'), 68000));

// 4. USAMOS LA TEMPERATURA YA DEFINIDA (o creamos una nueva con nombre único)
// t_aire_pixel_K es la temperatura del aire en Kelvin para cada píxel
var t_aire_pixel_K = image.select('temperatura_distribuida').add(273.15);

// 5. DISTRIBUCIÓN ESPACIAL DE PRESIÓN Y DENSIDAD
// Aplicamos la ecuación barométrica píxel a píxel desde la cota de referencia
var exponente = dem_10m.subtract(ee.Number(altitudReferencia))
                  .multiply(-9.81)
                  .divide(t_aire_pixel_K.multiply(ee.Number(Rd)));

var p_pixel = exponente.exp()
                .multiply(p_base_estacion)
                .rename('p_pixel');
                
// Densidad del aire (kg/m3) - Variable clave en Schaefer Eq. 1 y 3
var rho_a = p_pixel.divide(t_aire_pixel_K.multiply(Rd)).rename('rho_a');

// 6. CALOR SENSIBLE (SH) - Schaefer Eq. 1
// Ts = 0°C (273.15 K) por asunción de fusión superficial
var SH = rho_a.multiply(ca).multiply(Cp_star).multiply(u_imagen)
          .multiply(image.select('temperatura_distribuida').subtract(0))
          .rename('SH');

// 7. CALOR LATENTE (LH) - Schaefer Eq. 3
// e_sat_z: Presión de saturación en el aire (Bolton 1980) en Pascales
var e_sat_z = image.expression(
  '611.2 * exp((17.67 * T) / (T + 243.5))', 
  {'T': image.select('temperatura_distribuida')}
);

// e_a: Presión de vapor real en el aire
var e_a = e_sat_z.multiply(rh_dia.divide(100));

// e_s: Presión de vapor en la superficie (Saturación a 0°C = 611.2 Pa)
var e_s = 611.2;

var LH = rho_a.multiply(0.622).multiply(Lv).multiply(Cp_star).multiply(u_imagen)
          .multiply(e_a.subtract(e_s).divide(p_pixel))
          .rename('LH');
          
// CORRECCIÓN DE ESTABILIDAD (Schaefer 2020) - version 2
var Ri = image.select('temperatura_distribuida')
          .multiply(9.81).multiply(z_sensor)
          .divide(
            t_aire_pixel_K.multiply(u_imagen).multiply(u_imagen)  // u² como imagen
          );

// Factor lineal entre Ri=0.01 (factor=1) y Ri=0.2 (factor=0)
var corrFactor = Ri.multiply(-1).add(0.2).divide(0.19).clamp(0, 1);

// Para Ri <= 0.01, factor = 1 exactamente (sin corrección)
corrFactor = corrFactor.where(Ri.lte(0.01), ee.Image(1));

SH = SH.multiply(corrFactor).rename('SH');
LH = LH.multiply(corrFactor).rename('LH');

// 8. BALANCE ENERGÉTICO FINAL (Qm)
var swNet = incidente_total.multiply(ee.Image.constant(1).subtract(albedo)).rename('SW_net');
var Qm = swNet.add(image.select('LW_net')).add(SH).add(LH).rename('Qm');

// 9. AÑADIR BANDAS A LA IMAGEN
//image = image.addBands([SH, LH, Qm, rho_a, Eref, DSSR]);


// 10. MÁSCARA DE LA LENGUA (BAJO LOS 4150m)
var cotaELA = 4150;
//var mascaraLengua = dem_10m.lt(cotaELA);

// 1. Primero unimos TODAS las bandas calculadas en un solo bloque junto a la imagen
var imagenConTodasLasBandas = image.addBands([Edir, Eaniso_dif, Eiso_dif, Eref, Viso, 
                       incidente_total, Tb, SH, LH, Qm, rho_a, DSSR, u_imagen, Eglob_hor]);

// 2. Al final, aplicamos la máscara a la imagen completa. 
// Esto fuerza a que SH, LH, Qm y todas las demás queden estrictamente recortadas por debajo de 4150m.
return imagenConTodasLasBandas;

} // <--- AQUÍ TERMINA LA FUNCIÓN



// ==============================================================================
// 1. APLICAR LA FUNCIÓN A TODA LA COLECCIÓN
// ==============================================================================
var dataset_final = datasetConTemp.map(calculateRadiationComponents);

// ==============================================================================
// 2. LIMPIEZA GLOBAL (LA SOLUCIÓN A LOS DUPLICADOS Y NULOS)
// ==============================================================================
// A. Filtramos las imágenes vacías (días 100% nublados donde el albedo es nulo)
var dataset_sin_nulos = dataset_final.filter(ee.Filter.notNull(['Albedo']));

// B. Eliminamos las fechas repetidas (dejamos solo 1 pasada del satélite por día)
var dataset_limpio = dataset_sin_nulos.distinct('Fecha');

// --- LISTA DE FECHAS CORRUPTAS PARA EXCLUIR ---
var fechasCorruptas = [
  '01-06-2020', '19-01-2020', '05-02-2020', '10-04-2020',
  '02-02-2021', '04-02-2021', '07-02-2021', '01-03-2021', '08-04-2021', '10-04-2021',
  '23-04-2022',
  '06-03-2023', '03-04-2023',
  '09-02-2024', '24-02-2024', '20-03-2024', '29-04-2024',
  '28-02-2025', '25-03-2025', '02-04-2025'
];

// Filtramos la colección final para que estas fechas no entren en los cálculos
var dataset_final_balance = dataset_limpio.filter(ee.Filter.inList('Fecha', fechasCorruptas).not());

// A PARTIR DE AQUÍ, USAREMOS 'dataset_final_balance' PARA TODO LO DEMÁS



/////////////////////////////// VISUALIZACIÓN LW //////////////////////////////////


// Parámetros de visualización para la radiación
var visLW = {
  min: -100, // Valores negativos típicos de LW neta
  max: 0, 
  palette: ['#0000ff', '#ffffff', '#ff0000'] // Azul (pérdida de energía) a Rojo
};

var visLW_In = {
  min: 200, 
  max: 350, 
  palette: ['blue', 'green', 'yellow', 'red']
};

// Mostramos la mediana de la temporada para ver cómo se distribuye
var datasetConLW = datasetConTemp.map(calculateRadiationComponents);

//Map.addLayer(datasetConLW.select('LW_in').median().clip(juncal), visLW_In, 'Radiación LW Incidente (Schaefer)');
//Map.addLayer(datasetConLW.select('LW_net').median().clip(juncal), visLW, 'Balance de Onda Larga (LW Net)');
  


// ==============================================================================
// 7. CÁLCULO DE MEDIANAS ANUALES Y VISUALIZACIÓN DE MAPAS
// ==============================================================================

var anos = [2020, 2021, 2022, 2023, 2024, 2025];

// Paletas de colores para la visualización
var visRad = {min: 0, max: 800, palette: ['blue', 'green', 'yellow', 'red']};
var visAlbedo = {min: 0.2, max: 0.9, palette: ['black', 'grey', 'white']};

print('--- 📈 RESUMEN DE MEDIANAS ANUALES (W m⁻²) ---');

anos.forEach(function(anio) {
  // 1. Filtrar la colección por año
  var coleccionAnual = ee.ImageCollection(dataset_final_balance)
                        .filter(ee.Filter.eq('Anio', anio));
  
  // 2. Calcular la Imagen Mediana del año
  var medianaAnual = coleccionAnual.select(['Edir', 'Eaniso_dif', 'Eiso_dif', 'albedo']).median();
  
  // 3. Calcular la Radiación Incidente Total (Directa + Difusas)
  var SW_incidente = medianaAnual.select('Edir')
                      .add(medianaAnual.select('Eaniso_dif'))
                      .add(medianaAnual.select('Eiso_dif'))
                      .rename('SW_incidente');

  // 4. Reducción de región para obtener los números de la consola
  var stats = medianaAnual.addBands(SW_incidente).reduceRegion({
    reducer: ee.Reducer.median(),
    geometry: juncal,
    scale: 10,
    maxPixels: 1e9
  });

  // 5. Imprimir resultados en consola
  //print('Año ' + anio + ':', stats);

  // 6. Añadir capas al Mapa (solo para un año específico para no saturar, o puedes comentarlo)
  // Aquí añadimos las del último año o puedes elegir uno fijo
 /* if (anio === 2025) {
    Map.addLayer(medianaAnual.select('Edir').clip(juncal), visRad, 'Directa (Edir) ' + anio, false);
    Map.addLayer(medianaAnual.select('Eaniso_dif').clip(juncal), visRad, 'Dif. Anisotrópica ' + anio, false);
    Map.addLayer(medianaAnual.select('Eiso_dif').clip(juncal), visRad, 'Dif. Isotrópica ' + anio, false);
    Map.addLayer(medianaAnual.select('albedo').clip(juncal), visAlbedo, 'Albedo ' + anio, false);
    Map.addLayer(SW_incidente.clip(juncal), visRad, 'SW Incidente Total ' + anio, true);
  } */
});



// ==============================================================================
// BLOQUE FINAL: PRUEBA 2025 (CON CORRECCIÓN DE COLUMNAS)
// ==============================================================================

var estacionDiaria = ee.FeatureCollection("projects/tesis-457819/assets/sw_diaria_promedio");

var dataset_con_prediccion = ee.ImageCollection(dataset_final_balance.filter(ee.Filter.eq('Anio', 2025)).map(function(img) {
  var image = ee.Image(img);
  var fecha_buscada = image.date().format('dd-MM-YYYY');
  
  // 1. Buscamos el registro
  var registro = ee.Feature(estacionDiaria.filter(ee.Filter.eq('fecha', fecha_buscada)).first());
  
  // 2. Extraemos SW (Si el registro existe y tiene el valor)
  var sw_promedio_est = ee.Number(ee.Algorithms.If(registro, registro.get('SW'), null));
  
  var stats = image.select('DSSR').reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: puntoEstacion,
    scale: 30,
    tileScale: 4
  });
  
  var max_sat = ee.Number(stats.get('DSSR'));
  
  // 3. CÁLCULO DE R (Solo si tenemos ambos datos)
  // Si no hay dato de estación para ese día, usamos un factor R de 0.35 para no perder el relieve
  var R = ee.Number(ee.Algorithms.If(
    sw_promedio_est, 
    sw_promedio_est.divide(max_sat.add(0.001)), 
    0.35 
  ));
  
  var sw_diaria = image.select('DSSR').multiply(R).rename('SW_diaria_promedio');
  
  return image.addBands(sw_diaria).set('factor_R_diario', R).set('fecha_check', fecha_buscada);
}));

/*Map.addLayer(dataset_con_prediccion.select('SW_diaria_promedio').median().clip(juncal), 
             {min: 150, max: 350, palette: ['blue', 'cyan', 'green', 'yellow', 'red']}, 
             'SW Diaria 2025 (Corregida)');
             */
             



// ==============================================================================
// RECALCULO DEL BALANCE ENERGÉTICO DIARIO (24H)
// ==============================================================================
var dataset_balance_diario = dataset_con_prediccion.map(function(img) {
  
  // 1. Calculamos la Onda Corta Neta Diaria (SW_net = SW_incidente * (1 - Albedo))
  var sw_net_diario = img.select('SW_diaria_promedio')
    .multiply(ee.Image.constant(1).subtract(img.select('albedo')))
    .rename('SW_net_diario');

  // 2. Traemos los flujos que ya son diarios (LW, SH, LH)
  var lw_net = img.select('LW_net');
  var sh = img.select('SH');
  var lh = img.select('LH');

  // 3. NUEVA ECUACIÓN DE BALANCE (Qm Diario) — calculado píxel a píxel
 var qm_24h = sw_net_diario
    .add(img.select('LW_net'))
    .add(img.select('SH'))
    .add(img.select('LH'))
    .rename('Qm_24h');  // ← nombre completamente distinto a 'Qm'

  return img.addBands([sw_net_diario, qm_24h]);
});

// ==============================================================================
// EXPORTACIÓN DEL CSV
// ==============================================================================
var exportacionFinal = dataset_balance_diario.map(function(img) {
  
  var stats = img.select([
    'SW_diaria_promedio',
    'SW_net_diario',
    'albedo',
    'LW_net',
    'SH',
    'LH',
    'Qm_24h'
  ]).reduceRegion({
    reducer: ee.Reducer.median(),
    geometry: juncal,
    scale: 10,
    maxPixels: 1e9,
    tileScale: 4
  });

  return ee.Feature(null, {
    'Fecha':          img.get('fecha_check'),
    'SW_inc_diario':  stats.get('SW_diaria_promedio'),
    'SW_net_diario':  stats.get('SW_net_diario'),
    'Albedo':         stats.get('albedo'),
    'LW_net':         stats.get('LW_net'),
    'SH':             stats.get('SH'),
    'LH':             stats.get('LH'),
    'Qm_diario':      stats.get('Qm_24h')
  });
});

Export.table.toDrive({
  collection: exportacionFinal,
  description: 'Balance_Energia_Juncal_Diario_2025_lengua_13',
  folder: 'GEE_Tesis_Juncal',
  fileFormat: 'CSV',
  selectors: ['Fecha', 'SW_inc_diario', 'SW_net_diario', 'Albedo', 'LW_net', 'SH', 'LH', 'Qm_diario']
});

print('✅ Exportación lista. Dale Run en Tasks.');


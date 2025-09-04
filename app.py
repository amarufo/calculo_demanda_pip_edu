# Proyecto Flask: proyección poblacional con 4 datos y 3 tipos de proyección


import pandas as pd
import numpy as np
from scipy.stats import gmean
from flask import Flask, render_template, request, session, redirect, url_for
import math
import os



app = Flask(__name__)
app.secret_key = os.environ.get("SECRET_KEY", "ac33fb5ng54jhhld7???n33")




# Paso 1: Datos generales
@app.route("/", methods=["GET", "POST"])
def paso1():
    error = None
    datos = session.get('datos', {
        'nombre_proyecto': 'Nombre del Proyecto',
        'nombre_colegio': 'Nombre del Colegio',
        'nivel': 'Secundaria',
        'radio_influencia': 3,
        'area_distrito': 77.7,
        'est_by_aula': 30,
        'anio_form': 2025,
        'cantidad_anios_matricula': 5,
        'anio_i': 2027,
        'anio_f': 2036,
        "anio_censo1": 2007,
        "anio_censo2": 2017,
        "turnos": 2,
        "pob_censo1": 0,
        "pob_censo2": 0,
    })
    if request.method == 'POST':
        for campo in datos.keys():
            if campo in request.form:
                datos[campo] = request.form.get(campo, "")
        if datos.get("nivel") == "Primaria":
            datos["edades"] = list(range(6, 12))
        else:
            datos["edades"] = list(range(12, 17))
        datos["radio_influencia"] = float(datos["radio_influencia"])
        datos["area_distrito"] = float(datos["area_distrito"])
        datos["est_by_aula"] = float(datos["est_by_aula"])
        datos["anio_form"] = int(datos["anio_form"])
        datos["cantidad_anios_matricula"] = int(datos["cantidad_anios_matricula"])
        datos["anio_i"] = int(datos["anio_i"])
        datos["anio_f"] = int(datos["anio_f"])
        datos["anio_censo1"] = int(datos["anio_censo1"])
        datos["anio_censo2"] = int(datos["anio_censo2"])
        datos["turnos"] = int(datos["turnos"])
        # Creando nuevas variables
        anios_hist = list(range(datos["anio_form"] - datos["cantidad_anios_matricula"], datos["anio_form"]))
        datos["anios_hist"] = [int(anio) for anio in anios_hist]
        session["datos"] = datos
    # GET: NO guardar en sesión aquí, solo mostrar el formulario
    return render_template("paso1.html", datos=datos, error=error)

# Paso 2: Población, matrícula y no promovidos
@app.route("/paso2", methods=["GET", "POST"])
def paso2():
    error = None
    datos = session.get('datos', {})
    anio_form = datos.get("anio_form")
    cantidad_anios_matricula = datos.get("cantidad_anios_matricula")
    anio_censo1 = int(datos.get("anio_censo1"))
    anio_censo2 = int(datos.get("anio_censo2"))
    anios_hist = datos.get("anios_hist", [])
    edades = datos.get("edades", [])
    edades_ = [str(e) for e in edades]
    anios_proy = list(range(anio_form, anio_form + 5))
    anios_proy_ = [str(a) for a in anios_proy]
    anios_hist_ = [str(a) for a in anios_hist]

    if request.method == 'POST':
        # ---------------------------------------------------
        # post para guardar la población total por censos
        try:
            pob_censo1 = float(datos.get("pob_censo1", 0) or 0)
        except (TypeError, ValueError):
            pob_censo1 = 0.0
        try:
            pob_censo2 = float(datos.get("pob_censo2", 0) or 0)
        except (TypeError, ValueError):
            pob_censo2 = 0.0
        datos["pob_censo1"] = pob_censo1
        datos["pob_censo2"] = pob_censo2
        # ---------------------------------------------------
        # POST para guardar población por edades, matrícula y no promovidos
        dic_pop_edad_ = {}
        for anio in [anio_censo1, anio_censo2]:
            dic_pop_edad_[anio] = {}
            for edad in edades:
                try:
                    valor = int(request.form.get(f'pop_byedad{anio}_{edad}', 0) or 0)
                except (TypeError, ValueError):
                    valor = 0
                datos[f"dic_pop_edad_{anio}_{edad}"] = valor
                dic_pop_edad_[int(anio)][int(edad)] = valor
        datos["dic_pop_edad"] = dic_pop_edad_
        # ---------------------------------------------------
        dic_mat_by_anio_ = {}
        for anio in anios_hist:
            dic_mat_by_anio_[anio] = {}
            for edad in edades:
                try:
                    valor = int(request.form.get(f'matricula_{anio}_{edad}', 0) or 0)
                except (TypeError, ValueError):
                    valor = 0
                datos[f"dic_mat_by_anio_{anio}_{edad}"] = valor
                dic_mat_by_anio_[anio][edad] = valor
        datos["dic_mat_by_anio"] = dic_mat_by_anio_
        # ---------------------------------------------------
        dic_no_promv_ = {}
        for anio in anios_hist:
            dic_no_promv_[anio] = {}
            for edad in edades:
                try:
                    valor = int(request.form.get(f'noprom_{anio}_{edad}', 0) or 0)
                except (TypeError, ValueError):
                    valor = 0
                datos[f"dic_no_promv_{anio}_{edad}"] = valor
                dic_no_promv_[int(anio)][int(edad)] = valor
        datos["dic_no_promv"] = dic_no_promv_
        session["datos"] = datos
    # Guardar el diccionario correcto
    dic_ = {}
    for anio in [anio_censo1, anio_censo2]:
        dic_[anio] = {}
        for edad in edades:
            dic_[anio][edad] = datos.get(f"dic_pop_edad_{anio}_{edad}", 0)
    datos["dic_pop_edad"] = dic_
    dic_ = {}
    for anio in anios_hist:
        dic_[anio] = {}
        for edad in edades:
            dic_[anio][edad] = datos.get(f"dic_mat_by_anio_{anio}_{edad}", 0)
    datos["dic_mat_by_anio"] = dic_
    dic_ = {}
    for anio in anios_hist:
        dic_[anio] = {}
        for edad in edades:
            dic_[anio][edad] = datos.get(f"dic_no_promv_{anio}_{edad}", 0)
    datos["dic_no_promv"] = dic_
    session["datos"] = datos
    return render_template(
        "paso2.html",
        datos=datos,
        error=error
    )
# Paso 3: Proyección y gráfico
@app.route("/paso3", methods=["GET", "POST"])
def paso3():
    datos = session.get("datos", {})
    anio_form = datos.get("anio_form")
    cnt_anio_mat = datos.get("cantidad_anios_matricula")
    anio_f = datos.get("anio_f")
    anio_censo1 = int(datos.get("anio_censo1"))
    anio_censo2 = int(datos.get("anio_censo2"))
    pob_censo1 = float(datos.get("pob_censo1", 0)) if datos.get("pob_censo1") is not None else 0
    pob_censo2 = float(datos.get("pob_censo2", 0)) if datos.get("pob_censo2") is not None else 0
    anios_proy = list(range(anio_form - cnt_anio_mat, anio_f + 1))
    anios_hist = list(range(anio_form - cnt_anio_mat, anio_form))
    anio_i = int(datos.get("anio_i"))
    anio_f = int(datos.get("anio_f"))
    r = float(datos.get("radio_influencia"))
    area = float(datos.get("area_distrito", 0)) or 1  # Evitar división por cero
    dic_mat_by_anio = datos.get("dic_mat_by_anio", {})
    dic_no_promv = datos.get("dic_no_promv", {})
    edades = datos.get("edades", [])
    ###################################
    # CÁLCULO DE LA POBLACIÓN TOTAL
    ###################################
    dic_pob_total_ = {}
    try:
        tasa_poptotal = (pob_censo2 / pob_censo1) ** (1 / (anio_censo2 - anio_censo1)) - 1 
    except (ZeroDivisionError, TypeError):
        tasa_poptotal = 0
    for anio in anios_proy:
        if pob_censo2:
            dic_pob_total_[anio] = int(pob_censo2 * (1 + tasa_poptotal) ** (anio - anio_censo2))
        elif pob_censo1:
            dic_pob_total_[anio] = int(pob_censo1 * (1 + tasa_poptotal) ** (anio - anio_censo1))
        else:
            dic_pob_total_[anio] = 0
    datos["tasa_poptotal"] = tasa_poptotal
    datos["dic_pop_total"] = dic_pob_total_
    ###################################
    # CÁLCULO DE LA POBLACIÓN REFERENCIAL
    ###################################
    dic_pob_ref_ = {}
    for anio in anios_proy:
        dic_pob_ref_[anio] = int((math.pi * (r ** 2) / area) * dic_pob_total_[anio])
    datos["dic_pob_ref"] = dic_pob_ref_

    ###################################
    # CÁLCULO DE LA POBLACIÓN POTENCIAL
    ##################################
    anio_censo2_ = str(anio_censo2)
    anio_censo1_ = str(anio_censo1)
    edades_ = [str(e) for e in edades]
    dic_pop_edad = datos.get("dic_pop_edad", {})
    anios_proy_ = [str(a) for a in anios_proy]
    anios_hist_ = [str(a) for a in anios_hist]

    tasa_by_edad_ = {}
    for edad in edades_:
        a = dic_pop_edad[anio_censo1_][edad]
        b = dic_pop_edad[anio_censo2_][edad]
        if a > 0 and b > 0:
            tasa = (b / a) ** (1 / (anio_censo2 - anio_censo1)) - 1
        else:
            tasa = 0
        tasa_by_edad_[edad] = tasa
    datos["tasa_by_edad"] = tasa_by_edad_

    dic_pop_potencial = {}
    for edad in edades_:
        dic_pop_potencial[edad] = {}
        val_censo1 = dic_pop_edad[anio_censo1_][edad]
        val_censo2 = dic_pop_edad[anio_censo2_][edad]
        if val_censo2 > 0:
            base = val_censo2
            anio_base = anio_censo2
        else:
            base = val_censo1
            anio_base = anio_censo1
        pob_base = base * (math.pi * r**2) / area
        for anio, anio_ in zip(anios_proy, anios_proy_):
            dic_pop_potencial[edad][anio_] = int(pob_base * (1 + tasa_by_edad_[edad]) ** (anio - anio_base))
    datos["dic_pop_potencial"] = dic_pop_potencial
    ###################################
    # CÁLCULO DE LA EFECTIVA SIN PROYECTO
    ###################################
    # --- Tasa de proporción
    tasa_proporcion = {}
    for e in edades_:
        pop_edad = [dic_pop_potencial[e][a] for a in anios_hist_]
        mat_edad = [dic_mat_by_anio[a][e] for a in anios_hist_]
        prop = [m / p if p > 0 else 0 for p, m in zip(pop_edad, mat_edad)]
        tasa_proporcion[e] = gmean(prop)
    datos["tasa_proporcion"] = tasa_proporcion
    # --- Tasa de transición
    tasa_transicion = {}
    for e in edades[1:]:  # Desde la segunda edad hasta la última
        mat_up = [dic_mat_by_anio[str(a)][str(e)] for a in anios_hist[1:]]
        mat_down = [dic_mat_by_anio[str(a-1)][str(e-1)] for a in anios_hist[1:]]
        prop = [m / p if p > 0 else 0 for p, m in zip(mat_down, mat_up)]
        tasa_transicion[str(e)] = gmean(prop)
    datos["tasa_transicion"] = tasa_transicion
    # --- Proyección
    anio_proy = list(range(anio_form, anio_f + 1))
    anio_proy_ = [str(a) for a in anio_proy]


    e_min = min(edades)
    a_min = min(anio_proy)
    # Condiciones base
    pop_efec_sp = {}
    pop_efec_sp[str(e_min)] = {}
    for a in anio_proy_:
        dic = {a: dic_pop_potencial[str(e_min)][a] * tasa_proporcion[str(e_min)]}
        pop_efec_sp[str(e_min)].update(dic)

    edades2 = list(edades)
    edades2.remove(e_min)

    for e in edades2:
        pop_efec_sp[str(e)] = {}
        dic = {str(a_min): dic_pop_potencial[str(e)][str(a_min)] * tasa_proporcion[str(e)]}
        pop_efec_sp[str(e)].update(dic)
    # Completando la lista
    anios2 = list(anio_proy)
    anios2.remove(a_min)
    for e in edades2:
        for a in anios2:
            pop_efec_sp[str(e)][str(a)] = pop_efec_sp[str(e-1)][str(a-1)] * tasa_transicion[str(e)]
    # redondeando los valores
    for a in anio_proy:
        for e in edades:
            if pop_efec_sp[str(e)].get(str(a)) is None:
                pop_efec_sp[str(e)][str(a)] = 0
            else:
                pop_efec_sp[str(e)][str(a)] = pop_efec_sp[str(e)][str(a)]
    for a in anio_proy:
        for e in edades:
            try:
                valor = pop_efec_sp[str(e)][str(a)]
                pop_efec_sp[str(e)][str(a)] = round(float(valor) if valor not in [None, ""] else 0)
            except (TypeError, ValueError):
                pop_efec_sp[str(e)][str(a)] = 0
    datos["pop_efec_sp"] = pop_efec_sp
    ###################################
    # CÁLCULO DE LA EFECTIVA CON PROYECTO
    ###################################
    # Cálculo de parámetros
    proporcion_no_promovidos = {}
    for e in edades:
        no_prom = [dic_no_promv[str(a)][str(e)] for a in anios_hist]
        mat = [dic_mat_by_anio[str(a)][str(e)] for a in anios_hist]
        prop = [n / m if m > 0 else 0 for n, m in zip(no_prom, mat)]
        prop = [x for x in prop if x > 0]
        proporcion_no_promovidos[str(e)] = gmean(prop)
    datos["proporcion_no_promovidos"] = proporcion_no_promovidos
    tasas_cp = {}
    for e in edades:
        tasas_cp[str(e)] = tasa_proporcion[str(e)] * (1 + proporcion_no_promovidos[str(e)])
    datos["tasas_cp"] = tasas_cp
    # Crear el diccionario de población con proyecto
    anios_constr = list(range(anio_form, anio_i))
    anios_func = list(range(anio_i, anio_f+1))
    e_min = min(edades)
    edades3 = edades.copy()
    edades3.remove(e_min)
    pop_efec_cp = {}
    # condiciones base
    for e in edades:
        pop_efec_cp[str(e)] = {}
        for a in anio_proy:
            if a in anios_constr:
                pop_efec_cp[str(e)][str(a)] = pop_efec_sp[str(e)][str(a)]
            elif a in anios_func:
                if e == e_min:
                    pop_efec_cp[str(e)][str(a)] = dic_pop_potencial[str(e)][str(a)] * tasas_cp[str(e)]
                else:
                    try:
                        pop_efec_cp[str(e)][str(a)] = pop_efec_cp[str(e-1)][str(a-1)] * tasas_cp[str(e)]
                    except KeyError:
                        pop_efec_cp[str(e)][str(a)] = 0
    for a in anio_proy:
        for e in edades:
            if pop_efec_cp[str(e)].get(str(a)) is None:
                pop_efec_cp[str(e)][str(a)] = 0
            else:
                pass
    try:
        valor = pop_efec_cp[str(e)][str(a)]
        pop_efec_cp[str(e)][str(a)] = round(float(valor) if valor not in [None, ""] else 0)
    except (TypeError, ValueError):
        pop_efec_cp[str(e)][str(a)] = 0    

    datos["pop_efec_cp"] = pop_efec_cp

    ##################################
    # DIMENSIONAMIENTO - AULAS
    ##################################
    turnos = int(datos.get("turnos", 2))
    max_aulas_by_edad = {}
    aulas_necesarias_by_turno = {}

    for e in edades:
        # Calcula el máximo de aulas necesarias para esa edad
        lista = []
        for a in anio_proy:
            try:
                valor = pop_efec_cp[str(e)][str(a)] / (datos.get("est_by_aula", 30) if datos.get("est_by_aula", 30) > 0 else 30)
            except (TypeError, ValueError):
                valor = 0
            lista.append(valor) # Redondea hacia arriba
        max_aulas = max(lista)
        max_aulas_by_edad[str(e)] = max_aulas

        # Distribuye las aulas entre los turnos
        dic_aulas = {}
        base = max_aulas // turnos if turnos > 0 else max_aulas
        resto = max_aulas % turnos if turnos > 0 else 0
        for t in range(1, turnos + 1):
            dic_aulas[f'aulas_t{t}'] = base + (1 if t <= resto else 0)
        aulas_necesarias_by_turno[str(e)] = dic_aulas

    # Si quieres el máximo por turno para cada edad:
    aulas_necesarias_by_edad = {}
    for e in edades:
        aulas_por_turno = aulas_necesarias_by_turno[str(e)]
        max_aulas_turno = max(aulas_por_turno.values()) if aulas_por_turno else 0
        aulas_necesarias_by_edad[str(e)] = max_aulas_turno

    datos["aulas_necesarias_by_edad"] = aulas_necesarias_by_edad
    datos["aulas_necesarias_by_turno"] = aulas_necesarias_by_turno

                


    ##################################
    # ADAPTACIONES PARA EL RENDER
    ##################################
    # Tasa de crecimiento poblacional del distrito
    tasa_poptotal = datos.get("tasa_poptotal", 0)
    tic_dist_by_dep = {
        'AMAZONAS': {'tic_dist': 0.0024458725302197804},
        'ANCASH': {'tic_dist': -0.009722077727700852},
        'APURIMAC': {'tic_dist': -0.010176041095758207},
        'AREQUIPA': {'tic_dist': 0.002698241495252799},
        'AYACUCHO': {'tic_dist': -0.01915314237340786},
        'CAJAMARCA': {'tic_dist': -0.007530965006825317},
        'CALLAO': {'tic_dist': 0.007068794294313888},
        'CUSCO': {'tic_dist': -6.320814722067003e-05},
        'HUANCAVELICA': {'tic_dist': -0.027428666320475928},
        'HUANUCO': {'tic_dist': -0.02378710461378163},
        'ICA': {'tic_dist': 0.021343104369419066},
        'JUNIN': {'tic_dist': -0.005575224112039148},
        'LA LIBERTAD': {'tic_dist': 0.0020484522421301823},
        'LAMBAYEQUE': {'tic_dist': 0.01155271401947159},
        'LIMA': {'tic_dist': -0.000993094165232594},
        'LORETO': {'tic_dist': 0.001119823964029474},
        'MADRE DE DIOS': {'tic_dist': 0.03184191249719129},
        'MOQUEGUA': {'tic_dist': -0.020255832081772174},
        'PASCO': {'tic_dist': -0.013286263519479721},
        'PIURA': {'tic_dist': 0.0050009478848328705},
        'PUNO': {'tic_dist': -0.015083404010697627},
        'SAN MARTIN': {'tic_dist': 0.012149903727460838},
        'TACNA': {'tic_dist': -0.006558413043240586},
        'TUMBES': {'tic_dist': 0.02346913477418145},
        'UCAYALI': {'tic_dist': 0.0171539403116399},
        'TU PROYECTO': {'tic_dist': tasa_poptotal}
    }
    datos["tic_dist_by_dep"] = tic_dist_by_dep

    datos["area_influencia"] = round(math.pi * r**2, 2)
    turnos = int(datos.get("turnos", 2))
    list_turnos = list(range(1, turnos + 1))
    datos["list_turnos"] = list_turnos


    session["datos"] = datos
    return render_template(
        "paso3.html",
        datos=datos,
    )
if __name__ == "__main__":
    port = int(os.environ.get("PORT", 5000))
    app.run(debug=True, host="0.0.0.0", port=port)    
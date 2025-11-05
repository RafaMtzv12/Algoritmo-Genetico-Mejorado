# Algoritmo-Genetico-Mejorado
# 🧬 Algoritmo Genético Mejorado

**Autor:** Rafael de Jesus Martínez Velez 
**Matricula** 22056002
**Fecha:** Noviembre 2025  

---

## 🎯 Descripción
Este repositorio contiene una versión **propia y mejorada** del algoritmo genético para alineamiento de secuencias.  
La mejora se demuestra con un incremento del *fitness* respecto al algoritmo original, manteniendo la **validación de integridad**.

---

## ⚙️ Mejoras aplicadas
1. **Selección por torneo** → los mejores individuos tienen más probabilidad de reproducirse.  
2. **Elitismo** → conserva los mejores individuos de cada generación.  
3. **Mutación adaptativa** → ajusta automáticamente las probabilidades de mutación.  
4. **Mutación “gap-shift”** → mueve los gaps una posición para mejorar la alineación.  
5. **Validación de integridad** → garantiza que las secuencias sin gaps sigan siendo iguales.

---

## 📈 Evidencia
- `fitness_comparacion.png`: gráfica de comparación **original vs mejorado**  
- `fitness_comparacion.csv`: datos del *fitness* por generación  

---

## 🧪 Ejecución
```bash
pip install matplotlib pandas
python Algoritmo_mejorado.py
python compare_baseline_vs_mejora.py

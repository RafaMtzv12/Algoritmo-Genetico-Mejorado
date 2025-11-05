
# Runner para comparar el algoritmo original (AG10) y el mejorado
from Algoritmo_mejorado import correr as correr_mejorado
import random, copy, matplotlib.pyplot as plt, pandas as pd

# Para el baseline replicamos la lógica de AG10 (sin mutación por generación)
from Algoritmo_mejorado import (crear_poblacion_inicial, igualar_longitud, evaluar, 
                            cruzar_doble_punto_ignorar_gaps, eliminar_peores, get_sequences)

def correr_baseline(seed=7, gens=60, pop=10):
    random.seed(seed)
    pob=crear_poblacion_inicial(pop)
    # mutación inicial para introducir gaps
    from Algoritmo_mejorado import mutar_insertar_gaps
    pob=[mutar_insertar_gaps(ind,1,0.8) for ind in pob]
    pob=[igualar_longitud(ind) for ind in pob]
    scores=[evaluar(ind) for ind in pob]
    pob,scores=eliminar_peores(pob,scores)
    best_hist=[]
    originals=get_sequences()
    for g in range(gens):
        nueva=[]
        idx=list(range(len(pob))); random.shuffle(idx)
        if len(idx)%2==1: idx.append(idx[0])
        for i in range(0,len(idx),2):
            h1,h2=cruzar_doble_punto_ignorar_gaps(pob[idx[i]], pob[idx[i+1]])
            nueva.append(copy.deepcopy(pob[idx[i]]))
            nueva.append(copy.deepcopy(pob[idx[i+1]]))
            nueva.append(h1); nueva.append(h2)
        pob=nueva[:2*len(pob)]
        pob=[igualar_longitud(ind) for ind in pob]
        scores=[evaluar(ind) for ind in pob]
        pob,scores=eliminar_peores(pob,scores)
        assert all([ [a for a in s if a!='-']==[a for a in so if a!='-'] 
                     for ind in pob for s,so in zip(ind, originals)])
        best_hist.append(max(scores))
    return best_hist

def main():
    seed=7; gens=60
    hist_base=correr_baseline(seed, gens, pop=10)
    hist_mej, _, _ = correr_mejorado(seed, gens, pop=20, elitismo=2, torneo_k=3)
    df=pd.DataFrame({
        "generacion": list(range(1, gens+1)),
        "baseline_best": hist_base,
        "mejorado_best": hist_mej
    })
    df.to_csv("fitness_comparacion.csv", index=False)
    plt.figure()
    plt.plot(df["generacion"], df["baseline_best"], label="Original (AG10)")
    plt.plot(df["generacion"], df["mejorado_best"], label="Mejorado (elitismo+torneo+mutación adaptativa)")
    plt.xlabel("Generación"); plt.ylabel("Mejor fitness"); plt.title("Algoritmo original vs mejorado")
    plt.legend(); plt.tight_layout(); plt.savefig("fitness_comparacion.png")
    print("Guardado fitness_comparacion.csv y fitness_comparacion.png")

if __name__=="__main__":
    main()

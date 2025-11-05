
import random, copy, time
# BLOSUM62 table minimal (embedded)
BLOSUM62_STR = """
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V
A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0
R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3
N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3
D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3
C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1
Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2
E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2
G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3
H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3
I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3
L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1
K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2
M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1
F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1
P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2
S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0
W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3
Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1
V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4
"""
def parse_blosum62(s):
    rows=[r for r in s.strip().splitlines() if r.strip()]
    headers=rows[0].split()
    mat={}
    for r in rows[1:]:
        parts=r.split(); aa=parts[0]
        mat[aa]={h:int(v) for h,v in zip(headers, parts[1:])}
    return mat
B62=parse_blosum62(BLOSUM62_STR)

def get_sequences():
    seq1="MGSSHHHHHHSSGLVPRGSHMASMTGGQQMGRDLYDDDDKDRWGKLVVLGAVTQGQKLVVLGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQV"
    seq2="MKTLLVAAAVVAGGQGQAEKLVKQLEQKAKELQKQLEQKAKELQKQLEQKAKELQKQLEQKAKELQKQLEQKAGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQKELQKQLGQKAKEL"
    seq3="MAVTQGQKLVVLGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFAVVAGGQGQAEKLVKQLEQKAKELQKQLEQKAKELQKQLEQKAKELQKQLEQKAKELQKQLEQKALCVFAIN"
    return [list(seq1),list(seq2),list(seq3)]

def crear_poblacion_inicial(n=20):
    base=get_sequences()
    return [[row[:] for row in base] for _ in range(n)]

def igualar_longitud(individuo,gap='-'):
    m=max(len(f) for f in individuo)
    return [f+[gap]*(m-len(f)) for f in individuo]

def evaluar(individuo,gap_penal=-4):
    score=0; n=len(individuo); L=len(individuo[0])
    for c in range(L):
        for i in range(n):
            for j in range(i+1,n):
                a=individuo[i][c]; b=individuo[j][c]
                if a=='-' or b=='-': score+=gap_penal
                else: score+=B62[a][b]
    return score

def validar_poblacion(poblacion, originales):
    for ind in poblacion:
        for s,so in zip(ind,originales):
            if [a for a in s if a!='-'] != [a for a in so if a!='-']:
                return False
    return True

def validar_poblacion_detalle(poblacion, originales):
    ok = True
    for i, ind in enumerate(poblacion):
        for r, (s, so) in enumerate(zip(ind, originales)):
            a = [c for c in s if c != '-']
            b = [c for c in so if c != '-']
            if a != b:
                ok = False
                # busca primera diferencia
                m = min(len(a), len(b))
                mismatch_at = next((t for t in range(m) if a[t] != b[t]), None)
                print(f"[INTEGRIDAD] Falla en individuo {i}, fila {r}")
                print(f"  len(actual)={len(a)}, len(original)={len(b)}")
                if mismatch_at is not None:
                    t = mismatch_at
                    w0 = max(0, t-5); w1 = min(m, t+6)
                    print(f"  primer mismatch en pos {t}:")
                    print("  actual  :", ''.join(a[w0:w1]))
                    print("  original:", ''.join(b[w0:w1]))
                else:
                    print("  (mismo prefijo; longitudes distintas)")
    return ok


def mutar_insertar_gaps(individuo, n_gaps=1, p=0.3):
    nuevo=[]
    for fila in individuo:
        sec=fila[:]
        if random.random()<p:
            posiciones=set()
            for _ in range(n_gaps):
                pos=random.randint(0,len(sec))
                while pos in posiciones:
                    pos=random.randint(0,len(sec))
                posiciones.add(pos)
                sec.insert(pos,'-')
        nuevo.append(sec)
    return nuevo

def mutar_gap_shift(individuo, p=0.2):
    nuevo=[]
    for fila in individuo:
        sec=fila[:]
        if random.random()<p and '-' in sec:
            idxs=[i for i,a in enumerate(sec) if a=='-']
            i=random.choice(idxs)
            if i>0 and random.random()<0.5:
                sec[i],sec[i-1]=sec[i-1],sec[i]
            elif i<len(sec)-1:
                sec[i],sec[i+1]=sec[i+1],sec[i]
        nuevo.append(sec)
    return nuevo

def cruzar_doble_punto_ignorar_gaps(ind1, ind2):
    hijo1 = []
    hijo2 = []
    for s1, s2 in zip(ind1, ind2):
        # Quitamos gaps para trabajar solo con aminoácidos reales
        aaA = [a for a in s1 if a != '-']
        aaB = [a for a in s2 if a != '-']
        L = min(len(aaA), len(aaB))
        if L < 6:
            hijo1.append(s1[:])
            hijo2.append(s2[:])
            continue

        # Elegimos los puntos de corte sobre aminoácidos reales
        p1, p2 = sorted(random.sample(range(L), 2))

        # Construimos nuevas secuencias reales (sin gaps)
        nueva1 = aaA[:p1] + aaB[p1:p2] + aaA[p2:]
        nueva2 = aaB[:p1] + aaA[p1:p2] + aaB[p2:]

        # Reconstruimos con la estructura original de gaps
        def reconstruir(base, nueva):
            resultado = []
            k = 0
            for c in base:
                if c == '-':
                    resultado.append('-')
                else:
                    if k < len(nueva):
                        resultado.append(nueva[k])
                    else:
                        resultado.append('-')
                    k += 1
            return resultado

        hijo1.append(reconstruir(s1, nueva1))
        hijo2.append(reconstruir(s2, nueva2))
    return hijo1, hijo2

def eliminar_peores(poblacion, scores, porcentaje=0.5):
    idx=sorted(range(len(scores)), key=lambda i:scores[i], reverse=True)
    nsel=int(len(poblacion)*porcentaje)
    return [poblacion[i] for i in idx[:nsel]], [scores[i] for i in idx[:nsel]]

def torneo(poblacion, scores, k=3):
    idx=random.sample(range(len(poblacion)), k)
    best=max(idx, key=lambda i:scores[i])
    return copy.deepcopy(poblacion[best])

def correr(seed=7,gens=60,pop=20,elitismo=2,torneo_k=3):
    random.seed(seed)
    originals=get_sequences()
    # init
    pob=[mutar_insertar_gaps(ind,1,0.8) for ind in crear_poblacion_inicial(pop)]
    pob=[igualar_longitud(ind) for ind in pob]
    scores=[evaluar(ind) for ind in pob]
    mejor_actual=max(scores); estancado=0
    best_hist=[]
    for g in range(gens):
        if max(scores) <= mejor_actual: estancado+=1
        else: mejor_actual=max(scores); estancado=0
        p_ins=min(0.9, 0.3+0.05*estancado)
        p_shift=min(0.6, 0.1+0.03*estancado)
        # elitismo
        orden=sorted(range(len(scores)), key=lambda i:scores[i], reverse=True)
        elites=[copy.deepcopy(pob[i]) for i in orden[:elitismo]]
        # reproducción por torneo
        nueva=[]
        while len(nueva)<pop-elitismo:
            p1=torneo(pob, scores, torneo_k)
            p2=torneo(pob, scores, torneo_k)
            h1,h2=cruzar_doble_punto_ignorar_gaps(p1,p2)
            h1=mutar_insertar_gaps(h1,1,p_ins)
            h2=mutar_insertar_gaps(h2,1,p_ins)
            h1=mutar_gap_shift(h1,p_shift)
            h2=mutar_gap_shift(h2,p_shift)
            nueva.append(h1)
            if len(nueva)<pop-elitismo: nueva.append(h2)
        pob=elites+nueva[:pop-elitismo]
        pob=[igualar_longitud(ind) for ind in pob]
        scores=[evaluar(ind) for ind in pob]
        best_hist.append(max(scores))
        assert validar_poblacion(pob, originals)
    return best_hist, pob, scores
 
if __name__=="__main__":
    best_hist, pob, scores = correr()
    print("Mejor fitness final:", max(scores))
    originales = get_sequences()
    print("Validación de integridad:", validar_poblacion(pob, originales))


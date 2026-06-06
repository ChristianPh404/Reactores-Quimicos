import os

base_dir = r"f:\universidad\3º\REACTORES\2025\RQ\Reactores-Quimicos\RQ-master"
rcfp_tex_path = os.path.join(base_dir, "Capitulos", "RCFP.tex")

scripts = [
    {
        "name": "RCFP-1a",
        "file": "RCFP-1a.sce",
        "text": r"""\qs{RCFP-1a}{El equilibrio químico:
$$ A + B \rightleftharpoons C $$
Factor preexponencial de Arrhenius = $1.75 \times 10^8$ L/(mol$\cdot$h)\\
Factor preexponencial de Van't Hoff = $8.25 \times 10^{-22}$ L/mol\\
Energía de activación = $62350$ J/mol\\
Entalpía de reacción = $-136400$ J/mol

Se lleva a cabo en condiciones isotermas a 310 K en un reactor tubular de 800 dm de longitud y 3 dm de diámetro alimentado por 50 L/h de los compuestos A, B y C a concentración 1.5, 2.0 y 0.1 mol/L, respectivamente. La mezcla reaccionante presenta una densidad de 1150 g/L y una capacidad calorífica de 3.8 J/(g$\cdot$K).

\begin{enumerate}[label={\alph*)}]
    \item Obtener el estado estacionario del proceso.
\end{enumerate}}
"""
    },
    {
        "name": "RCFP-1b",
        "file": "RCFP-1b.sce",
        "text": r"""\qs{RCFP-1b}{El equilibrio químico:
$$ A + B \rightleftharpoons C $$
Factor preexponencial de Arrhenius = $1.75 \times 10^8$ L/(mol$\cdot$h)\\
Factor preexponencial de Van't Hoff = $8.25 \times 10^{-22}$ L/mol\\
Energía de activación = $62350$ J/mol\\
Entalpía de reacción = $-136400$ J/mol

Se lleva a cabo en condiciones adiabáticas en un reactor tubular de 800 dm de longitud y 3 dm de diámetro alimentado por 50 L/h a 310 K de los compuestos A, B y C a concentración 1.5, 2.0 y 0.1 mol/L, respectivamente. La mezcla reaccionante presenta una densidad de 1150 g/L y una capacidad calorífica de 3.8 J/(g$\cdot$K).

\begin{enumerate}[label={\alph*)}]
    \item Obtener el estado estacionario del proceso.
\end{enumerate}}
"""
    },
    {
        "name": "RCFP-1c",
        "file": "RCFP-1c.sce",
        "text": r"""\qs{RCFP-1c}{El equilibrio químico:
$$ A + B \rightleftharpoons C $$
Factor preexponencial de Arrhenius = $1.75 \times 10^8$ L/(mol$\cdot$h)\\
Factor preexponencial de Van't Hoff = $8.25 \times 10^{-22}$ L/mol\\
Energía de activación = $62350$ J/mol\\
Entalpía de reacción = $-136400$ J/mol

Se lleva a cabo en condiciones no adiabáticas en un reactor tubular de 800 dm de longitud y 3 dm de diámetro alimentado por 50 L/h a 310 K de los compuestos A, B y C a concentración 1.5, 2.0 y 0.1 mol/L, respectivamente. El reactor está refrigerado por una camisa por la que circula un fluido a 310 K, con un coeficiente global de transmisión de calor de 900 J/(dm$^2\cdot$h$\cdot$K). La mezcla reaccionante presenta una densidad de 1150 g/L y una capacidad calorífica de 3.8 J/(g$\cdot$K).

\begin{enumerate}[label={\alph*)}]
    \item Obtener el estado estacionario del proceso.
\end{enumerate}}
"""
    },
    {
        "name": "RCFP-1d",
        "file": "RCFP-1d.sce",
        "text": r"""\qs{RCFP-1d}{El equilibrio químico:
$$ A + B \rightleftharpoons C $$
Factor preexponencial de Arrhenius = $1.75 \times 10^8$ L/(mol$\cdot$h)\\
Factor preexponencial de Van't Hoff = $8.25 \times 10^{-22}$ L/mol\\
Energía de activación = $62350$ J/mol\\
Entalpía de reacción = $-136400$ J/mol

Se lleva a cabo en un sistema de dos reactores tubulares adiabáticos en serie con un enfriamiento intermedio. El sistema completo tiene una longitud equivalente de 800 dm y 3 dm de diámetro, y es alimentado por 50 L/h a 310 K de los compuestos A, B y C a concentración 1.5, 2.0 y 0.1 mol/L, respectivamente. Entre ambos reactores se enfría la mezcla hasta la temperatura inicial de 310 K. La mezcla reaccionante presenta una densidad de 1150 g/L y una capacidad calorífica de 3.8 J/(g$\cdot$K).

\begin{enumerate}[label={\alph*)}]
    \item Obtener el estado estacionario del proceso y optimizar la longitud del primer reactor.
\end{enumerate}}
"""
    },
    {
        "name": "RCFP-2",
        "file": "RCFP-2.sce",
        "text": r"""\qs{RCFP-2}{La reacción elemental e isoterma:
$$ A + B \rightarrow 2B $$
Con constante cinética $k = 1$ L/(mol$\cdot$h), se lleva a cabo en un proceso con recirculación. El sistema consta de un reactor tubular de 12 dm de longitud y 1 dm de diámetro. La alimentación fresca es de 5 L/h con una concentración de 2 mol/L de A y 0.01 mol/L de B. Parte de la corriente de salida se recircula a la entrada del reactor, con una relación de recirculación de $R = 0.1$.

\begin{enumerate}[label={\alph*)}]
    \item Obtener el estado estacionario del proceso con el método del punto fijo.
    \item Optimizar la fracción de recirculación.
\end{enumerate}}
"""
    }
]

with open(rcfp_tex_path, "r", encoding="utf-8") as f:
    content = f.read()

# Replace \section{Scilab} with the content
target = r"\section{Scilab}"
replacement = target + "\n\n"

for s in scripts:
    with open(os.path.join(base_dir, "3-RCFP", s["file"]), "r", encoding="utf-8") as f:
        code = f.read()
    
    # ensure proper encoding of special characters is left intact
    replacement += s["text"]
    replacement += "\\begin{minted}[fontsize=\\small, linenos]{scilab}\n"
    replacement += code
    if not replacement.endswith("\n"):
        replacement += "\n"
    replacement += "\\end{minted}\n\n"

content = content.replace(target, replacement)

with open(rcfp_tex_path, "w", encoding="utf-8") as f:
    f.write(content)

print("Successfully updated RCFP.tex with all scripts!")

import os

base_dir = r"f:\universidad\3º\REACTORES\2025\RQ\Reactores-Quimicos\RQ-master"
rdmp_tex_path = os.path.join(base_dir, "Capitulos", "RDMP.tex")

table_code = r"""
\vspace{1em}
\begin{table}[H]
    \centering
    \renewcommand{\arraystretch}{1.5}
    \resizebox{\textwidth}{!}{
    \begin{tabular}{|p{4.5cm}|p{2.5cm}|p{2cm}|p{2.5cm}|p{2.5cm}|p{3.5cm}|}
        \hline
        \rowcolor{gray!20}
        \textbf{Objetivo / Detalle Clave} & \textbf{Isotermo} & \textbf{Adiabático} & \textbf{No adiabático (T camisa cte)} & \textbf{No adiabático (T camisa var)} & \textbf{Exámenes y Otros} \\
        \hline
        \textbf{Calcular conversión final / dinámica} & 
        \hyperref[qs:RDMP-1]{RDMP-1}\newline \hyperref[qs:RDMP-MULT-1]{RDMP-MULT-1} & 
        \hyperref[qs:RDMP-2]{RDMP-2}\newline \hyperref[qs:RDMP-5]{RDMP-5} & 
        \hyperref[qs:RDMP-2022]{RDMP-2022}\newline \hyperref[qs:RDMP-MULTI-2024-0]{RDMP-M-2024} & 
        & 
        \hyperref[qs:RDMP-6]{RDMP-6}\newline \hyperref[qs:RDMP-O-2018]{RDMP-O-2018}\newline \hyperref[qs:RDMP-O-2022]{RDMP-O-2022} \\
        \hline
        \textbf{Encontrar extremo (Tmax, Tmin, o máx. conc)} & 
        \hyperref[qs:RDMP-4]{RDMP-4}\newline \hyperref[qs:RDMP-MULT-1]{RDMP-MULT-1} & 
        \hyperref[qs:RDMP-2]{RDMP-2} & 
        \hyperref[qs:RDMP-MULTI-2024-0]{RDMP-M-2024}\newline \hyperref[qs:RDMP-MULTI-2025-O]{RDMP-M-2025} & 
        \hyperref[qs:RDMP-3b]{RDMP-3b} & 
        \hyperref[qs:RDMP-2023-O]{RDMP-2023-O}\newline \hyperref[qs:RDMP-3a]{RDMP-3a}\newline \hyperref[qs:RDMP-3c]{RDMP-3c}\newline \hyperref[qs:RDMP-6]{RDMP-6}\newline \hyperref[qs:RDMP-E-2020]{RDMP-E-2020}\newline \hyperref[qs:RDMP-E-2022]{RDMP-E-2022} \\
        \hline
        \textbf{Determinar diseño térmico (caudal/T agua)} & 
        & 
        & 
        \hyperref[qs:RDMP-2022]{RDMP-2022}\newline \hyperref[qs:RDMP-MULT-2]{RDMP-MULT-2}\newline \hyperref[qs:RDMP-MULTI-2024-0]{RDMP-M-2024}\newline \hyperref[qs:RDMP-MULTI-2025-O]{RDMP-M-2025} & 
        \hyperref[qs:RDMP-3b]{RDMP-3b} & 
        \hyperref[qs:RDMP-2023-O]{RDMP-2023-O}\newline \hyperref[qs:RDMP-3a]{RDMP-3a}\newline \hyperref[qs:RDMP-E-2020]{RDMP-E-2020}\newline \hyperref[qs:RDMP-E-2022]{RDMP-E-2022}\newline \hyperref[qs:RDMP-O-2018]{RDMP-O-2018}\newline \hyperref[qs:RDMP-O-2022]{RDMP-O-2022} \\
        \hline
        \textbf{Localizar donde se igualan concentraciones} & 
        & 
        & 
        & 
        & 
        \hyperref[qs:RDMP-O-2022]{RDMP-O-2022} \\
        \hline
        \textbf{Progresión óptima de temperatura} & 
        & 
        & 
        & 
        & 
        \hyperref[qs:RDMP-3c]{RDMP-3c} \\
        \hline
    \end{tabular}
    }
    \caption{Clasificación de ejercicios del RDMP por tipo de operación y objetivos clave.}
    \label{tab:clasificacion_rdmp}
\end{table}
\vspace{1em}
"""

with open(rdmp_tex_path, "r", encoding="utf-8") as f:
    content = f.read()

target = r"\section{Problemas Practicos-Scilab}"
replacement = target + "\n" + table_code

content = content.replace(target, replacement)

with open(rdmp_tex_path, "w", encoding="utf-8") as f:
    f.write(content)

print("Table successfully inserted!")

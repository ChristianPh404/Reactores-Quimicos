import os

base_dir = r"f:\universidad\3º\REACTORES\2025\RQ\Reactores-Quimicos\RQ-master"
rcmp_tex_path = os.path.join(base_dir, "Capitulos", "RCMP.tex")

table_code = r"""
\vspace{1em}
\begin{table}[H]
    \centering
    \renewcommand{\arraystretch}{1.5}
    \resizebox{\textwidth}{!}{
    \begin{tabular}{|p{4.5cm}|p{2.5cm}|p{3.0cm}|p{3.0cm}|p{3.5cm}|}
        \hline
        \rowcolor{gray!20}
        \textbf{Objetivo / Detalle Clave} & \textbf{Isotermo} & \textbf{No adiabático (T camisa cte)} & \textbf{No adiabático (T camisa var)} & \textbf{Exámenes y Otros} \\
        \hline
        \textbf{Estudiar dinámica / Conversión básica} & 
        \hyperlink{qs:RCMP-1}{RCMP-1}\newline \hyperlink{qs:RCMP-2}{RCMP-2} & 
        & 
        & 
        \\
        \hline
        \textbf{Optimización / Maximizar selectividad} & 
        \hyperlink{qs:RCMP-MULT-1}{RCMP-MULT-1}\newline \hyperlink{qs:RCMP-4}{RCMP-4} & 
        & 
        & 
        \\
        \hline
        \textbf{Diseño térmico (caudal/temperatura agua)} & 
        & 
        & 
        \hyperlink{qs:RCMP-3}{RCMP-3}\newline \hyperlink{qs:RCMP-MULT-2}{RCMP-MULT-2} & 
        \hyperlink{qs:RCMP-2023-E}{RCMP-2023-E}\newline \hyperlink{qs:RCMP-2023-O}{RCMP-2023-O}\newline \hyperlink{qs:RCMP-MULTI-2024-E}{RCMP-MULTI-2024-E}\newline \hyperlink{qs:RCMP-MULTI-2025-E}{RCMP-MULTI-2025-E} \\
        \hline
        \textbf{Estados estacionarios y Estabilidad} & 
        & 
        \hyperlink{qs:RCMP-5}{RCMP-5} & 
        \hyperlink{qs:RCMP-3}{RCMP-3}\newline \hyperlink{qs:RCMP-MULT-2}{RCMP-MULT-2} & 
        \hyperlink{qs:RCMP-2023-E}{RCMP-2023-E} \\
        \hline
    \end{tabular}
    }
    \caption{Clasificación de ejercicios del RCMP por tipo de operación y objetivos clave.}
    \label{tab:clasificacion_rcmp}
\end{table}
\vspace{1em}
"""

with open(rcmp_tex_path, "r", encoding="utf-8") as f:
    content = f.read()

target = r"\section{Problemas Practicos-Scilab}"
replacement = target + "\n" + table_code

if "\\label{tab:clasificacion_rcmp}" not in content:
    content = content.replace(target, replacement)

    with open(rcmp_tex_path, "w", encoding="utf-8") as f:
        f.write(content)

    print("Table successfully inserted into RCMP.tex!")
else:
    print("Table already exists in RCMP.tex!")

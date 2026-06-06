import os, glob, re

scripts = glob.glob('2-RCMP/*.sce')
for s in scripts:
    with open(s, encoding='utf-8') as f:
        content = f.read()
    
    op_type = 'Desconocido'
    if 'Isotermo' in content or 'isotermo' in content.lower(): op_type = 'Isotermo'
    elif 'Adiabático' in content or 'Adiabatico' in content or 'adiabático' in content.lower(): op_type = 'Adiabático'
    
    if 'No adiabático' in content or 'No adiabatico' in content or 'no adiabático' in content.lower():
        if 'camisa' in content.lower() and 'variable' in content.lower():
            op_type = 'No adiabático (T camisa var)'
        elif 'serpent' in content.lower():
            op_type = 'No adiabático (Serpentín)'
        else:
            op_type = 'No adiabático (T camisa cte)'
            
    details = []
    if 'conversion' in content.lower() or 'conversión' in content.lower() or 'xa' in content.lower(): details.append('Calcular conversión')
    if 'max(' in content or 'min(' in content: details.append('Encontrar extremo (Tmax/Tmin)')
    if 'optimo' in content.lower() or 'óptimo' in content.lower(): details.append('Progresión óptima')
    if 'caudal' in content.lower() or 'tj' in content.lower(): details.append('Determinar caudal/temperatura agua (Tj)')
    if 'igual' in content.lower(): details.append('Localizar donde se igualen concentraciones')
    if 'selectivi' in content.lower(): details.append('Calcular selectividad')
    if 'estabilidad' in content.lower(): details.append('Analizar estabilidad')
    
    print(f'{os.path.basename(s).replace(".sce", "")} | {op_type} | {", ".join(details)}')

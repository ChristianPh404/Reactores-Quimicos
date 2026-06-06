import re

def process_file(filepath):
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()

    # Find all \qs blocks and their positions
    # We will split the file by \qs{
    # But \qs{ can be anywhere. Let's use a regex to match the \qs{...} lines exactly.
    # Actually, a better way is to split by \qs{ and then find the matching closing brace for the second argument.
    
    parts = re.split(r'(\\qs\{[^}]+\}\{)', content)
    
    if len(parts) == 1:
        return # No \qs found
        
    # We will try to parse blocks by finding \qs{Name}{
    
    blocks = []
    
    lines = content.split('\n')
    
    prefix_lines = []
    qs_blocks = []
    
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith('\\qs{'):
            # start of a block
            name_match = re.match(r'\\qs\{([^}]+)\}\{', line)
            if not name_match:
                prefix_lines.append(line)
                i += 1
                continue
            
            name = name_match.group(1)
            
            # extract the entire block by counting braces
            block_lines = [line]
            open_braces = line.count('{') - line.count('}')
            
            i += 1
            while i < len(lines) and open_braces > 0:
                block_lines.append(lines[i])
                open_braces += lines[i].count('{') - lines[i].count('}')
                i += 1
            
            qs_blocks.append({
                'name': name,
                'lines': block_lines
            })
            
        else:
            if not qs_blocks:
                prefix_lines.append(line)
            else:
                # If we have already seen a qs block, but there is some text between qs blocks
                # We append it to the previous qs_block
                qs_blocks[-1]['lines'].append(line)
            i += 1

    # Now filter the qs_blocks that are "dated"
    # For RDMP, dated ones have years like 2020, 2021, or ?
    # Let's say any block with "-O-" or "-E-" or "20" in the name or "?" is considered a dated exam if it's towards the end.
    
    # We'll identify dated blocks
    dated_blocks = []
    undated_blocks = []
    
    for block in qs_blocks:
        name = block['name']
        
        # rename logic: RDMP-O-2023 -> RDMP-2023-O
        m = re.match(r'([A-Z]+)-(O|E)-(\d{4}|\?)', name)
        if m:
            prefix, term, year = m.groups()
            new_name = f"{prefix}-{year}-{term}"
            # Replace in the first line
            block['lines'][0] = block['lines'][0].replace(f"\\qs{{{name}}}{{", f"\\qs{{{new_name}}}{{")
            block['name'] = new_name
            block['year'] = year
            block['term'] = term
            dated_blocks.append(block)
        else:
            m2 = re.match(r'([A-Z]+)-(\?|2\d{3})-(E|O|\?)', name)
            if m2:
                prefix, year, term = m2.groups()
                block['year'] = year
                block['term'] = term
                dated_blocks.append(block)
            else:
                m3 = re.match(r'([A-Z]+)-\?-\?', name)
                if m3:
                    prefix = m3.group(1)
                    new_name = f"{prefix}-?-?"
                    block['lines'][0] = block['lines'][0].replace(f"\\qs{{{name}}}{{", f"\\qs{{{new_name}}}{{")
                    block['name'] = new_name
                    block['year'] = '?'
                    block['term'] = '?'
                    dated_blocks.append(block)
                else:
                    undated_blocks.append(block)

    if not dated_blocks:
        print(f"No dated blocks found in {filepath}")
        return

    # Sort dated blocks
    # sort key: year first (with ? = '0000'), then term
    def sort_key(b):
        y = b['year']
        if y == '?': y = '0000'
        t = b['term']
        if t == '?': t = 'A'
        return (y, t)

    dated_blocks.sort(key=sort_key)

    # Reconstruct the file
    with open(filepath, 'w', encoding='utf-8') as f:
        f.write('\n'.join(prefix_lines))
        f.write('\n')
        
        for b in undated_blocks:
            f.write('\n'.join(b['lines']))
            f.write('\n')
            
        for b in dated_blocks:
            f.write('\n'.join(b['lines']))
            f.write('\n')
            
    print(f"Processed {filepath} successfully.")

process_file('Capitulos/RDMP.tex')
process_file('Capitulos/RCFP.tex')
process_file('Capitulos/RSMP.tex')

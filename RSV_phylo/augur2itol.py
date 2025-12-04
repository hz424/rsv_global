import json
from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import to_rgba
from matplotlib.patches import Patch
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import random
from collections import Counter
import pdb

asMuts = [42, 205, 206, 209, 211, 220, 255, 312, 424, 467] ## loci in antigenic sites were mutations are observed

# see "RSV presentation_edits_hz_ah.pptx"
haoMuts = [12, 15, 19, 42, 103, 113, 116, 127, 173, 190, 191, 206, 209, 211, 220, 327, 389, 466, 467, 529] # RSV B 
haoMuts = [3, 12, 14, 15, 19, 20, 23, 103, 119, 122, 276, 381, 419, 518] # RSV A
entryThreshold = 0.9 ## threshold for origin probability 

natDict = {'EGYPTIAN': 'NORTH AFRICAN', 'SUDANESE': 'NORTH AFRICAN', 'MOROCCAN': 'NORTH AFRICAN', 'MAURITANIAN': 'NORTH AFRICAN', 
 'SOMALI': 'EAST AFRICAN', 'ETHIOPIAN': 'EAST AFRICAN', 'LEBANESE': 'LEVANT/TURKEY', 'PALESTINIAN': 'LEVANT/TURKEY', 
 'SYRIAN': 'LEVANT/TURKEY', 'IRAQI': 'LEVANT/TURKEY', 'TURKISH': 'LEVANT/TURKEY', 'INDIAN': 'SOUTH EAST ASIAN', 'PAKISTANI': 'SOUTH EAST ASIAN', 'BANGLADESHI': 
 'SOUTH EAST ASIAN', 'SRI LANKAN': 'SOUTH EAST ASIAN', 'YEMENI': 'ARABIA', 'OMANI': 'ARABIA', 'SAUDI': 'ARABIA', 'AMERICAN': 'NORTH AMERICAN', 'CANADIAN': 'NORTH AMERICAN', 'EMIRATI': 'EMIRATI', 
 'TANZANIAN': 'AFRICAN', 'UGANDAN': 'AFRICAN', 'SOUTH AFRICAN': 'AFRICAN', 'JORDANIAN': 'LEVANT/TURKEY'}

def get_hex_color_from_colormap(value, colormap_name='viridis', alpha=1.0):
    cmap = cm.get_cmap(colormap_name)
    rgba_color = cmap(value, alpha=alpha)
    return mcolors.to_hex(rgba_color, keep_alpha=(alpha < 1.0))



color_palette = {
#    "blue_shades": {
        "A.D.5.2": "#ADD8E6",
        "A.D.5.1": "#4682B4",
        "A.D.5.3": "#191970",
#    "red_shades": {
        "A.D.3": "#F08080",
        "A.D.3.1": "#DC143C",
        "A.D.3.3": "#8B0000",
#    "yellow_shades": {
        "A.D.1": "#FFFFED",
        "A.D.1.4": "#FFD700",
        "A": "#F0F0F0"
}

import colorsys

def generate_blue_shades(hue = 220, n=7):
    """
    Generate n maximally different shades of blue.
    Returns a list of RGB hex color codes.
    """
    hue /= 360
    # Vary saturation and value to maximize perceptual distance
    shades = []
    for i in range(n):
        # Alternate brightness and saturation to maximize distance
        saturation = 0.8 # + (i % 2) * 0.7  # alternate low/high saturation
        value = 0.4 + (i / (n - 1)) * 0.6  # evenly spread value from dark to bright

        # Convert HSV -> RGB (0..1 range)
        r, g, b = colorsys.hsv_to_rgb(hue, saturation, value)

        # Convert to hex
        hex_color = "#{:02x}{:02x}{:02x}".format(int(r * 255), int(g * 255), int(b * 255))
        shades.append(hex_color)

    return shades

def gray_scale(steps=10):
    import numpy as np
    return [f"#{int(a):02x}{int(a):02x}{int(a):02x}" for a in np.linspace(0, 255, steps)]

def draw_pie(ax, ratios, X, Y, size, colors=None, alpha=1.0):
    rgba_colors = [to_rgba(c, alpha=alpha) for c in colors[:len(ratios)]]
    wedges, _ = ax.pie(
        ratios,
        radius=size,
        colors=rgba_colors,
        center=(X, Y+size),
        wedgeprops={'edgecolor': 'black'},
        startangle=0
    )

def collapseMutHist(Fmuts):
    """removes (circular) backmutations"""
    from collections import defaultdict
    fd = defaultdict(list)
    for f in Fmuts:
        key = int(f[1:-1])
        fd[key].append((f[0],f[-1])) 
    return {pos:mhist[-1][-1] for pos, mhist in fd.items() if mhist[0][0]!=mhist[-1][-1]}   

def antigenicSiteMutScores(mutations, binary=True):
    mutScore = lambda sr: min(sum(1 for locus in sr if locus in mutations), 1)
    if not binary:
        mutScore = lambda sr: sum(1 for locus in sr if locus in mutations)/len(sr)
    return {site: mutScore(siteRange) for site, siteRange in sites.items()}

def calcAntigenicDrift(mutations):
    '''1. given mutations (list, diverg from ref): calculate full sequence, using ref (KT992...)
       2. calc. mutations wrt to Nirsevimab reference KX ...'''
    seq = list(ref)
    # {105: 'S', 574: 'N', 276: 'S', 540: 'A', 12: 'I', 169: 'N'} # 
    for key, val in mutations.items():
        if key <= len(seq):
            seq[key-1] = val
    return {i+1: res2 for i, (res1, res2) in enumerate(zip(nirRefA, seq)) if res1!=res2}
        
subcladeCounter = 0
verbose = True
subtype = 'A'
baseDir = "/home/ahenschel/Flu_Madikay/RSV"

nirRefA = """MELPILKTNAITTILAAVTLCFASSQNITEEFYQSTCSAVSKGY
LSALRTGWYTSVITIELSNIKENKCNGTDAKVKLIKQELDKYKNAVTELQLLMQSTPA
ANSRARRELPRFMNYTLNNTKNTNVTLSKKRKRRFLGFLLGVGSAIASGIAVSKVLHL
EGEVNKIKSALLSTNKAVVSLSNGVSVLTSKVLDLKNYIDKQLLPIVNKQSCSISNIE
TVIEFQQKNNRLLEITREFSVNAGVTTPVSTYMLTNSELLSLINDMPITNDQKKLMSS
NVQIVRQQSYSIMSIIKEEVLAYVVQLPLYGVIDTPCWKLHTSPLCTTNTKEGSNICL
TRTDRGWYCDNAGSVSFFPQAETCKVQSNRVFCDTMNSLTLPSEVNLCNIDIFNPKYD
CKIMTSKTDVSSSVITSLGAIVSCYGKTKCTASNKNRGIIKTFSNGCDYVSNKGVDTV
SVGNTLYYVNKQEGKSLYVKGEPIINFYDPLVFPSDEFDASISQVNEKINQSLAFIRK
SDELLHNVNAGKSTTNIMITTIIIVIIVILLALIAVGLLLYCKARSTPVTLSKDQLSG
INNIAFSN""".replace('\n','')

ref = '''MELLILKANAITTILTAVTFCFASGQNITEEFYQSTCSAVSKGYL
                     SALRTGWYTSVITIELSNIKKNKCNGTDAKVKLIKQELDKYKNAVTELQLLMQSTQATN
                     NRARRELPRFMNYTLNNAKKTNVTLSKKRKRRFLGFLLGVGSAIASGVAVSKVLHLEGE
                     VNKIKSALLSTNKAVVSLSNGVSVLTSKVLDLKNYIDKQLLPIVNKQSCSISNIETVIE
                     FQQKNNRLLEITREFSVNAGVTTPVSTYMLTNSELLSLINDMPITNDQKKLMSNNVQIV
                     RQQSYSIMSIIKEEVLAYVVQLPLYGVIDTPCWKLHTSPLCTTNTKEGSNICLTRTDRG
                     WYCDNAGSVSFFPQAETCKVQSNRVFCDTMNSLTLPSEVNLCNVDIFNPKYDCKIMTSK
                     TDVSSSVITSLGAIVSCYGKTKCTASNKNRGIIKTFSNGCDYVSNKGVDTVSVGNTLYY
                     VNKQEGKSLYVKGEPIINFYDPLVFPSDEFDASISQVNEKINQSLAFIRKSDELLHNVN
                     AGKSTTNIMITTIIIVIIVILLSLIAVGLLLYCKARSTPVTLSKDQLSGINNIAFSN'''.replace('\n','').replace(' ','')

if subtype == 'B':
    nirRefA = '''MELLIHRSSAIFLTFAINALYLTSSQNITEEFYQSTCSAVSRGY
                        LSALRTGWYTSVITIELSNIKETKCNGTDTKVKLIKQELDKYKNAVTELQLLMQNTPA
                        ANNRARREAPQYMNYTINTTKNLNVSISKKRKRRFLGFLLGVGSAIASGIAVSKVLHL
                        EGEVNKIKNALLSTNKAVVSLSNGVSVLTSKVLDLKNYINNQLLPIVNQQSCRISNIE
                        TVIEFQQKNSRLLEITREFSVNAGVTTPLSTYMLTNSELLSLINDMPITNDQKKLMSS
                        NVQIVRQQSYSIMSIIKEEVLAYVVQLPIYGVIDTPCWKLHTSPLCTTNIKEGSNICL
                        TRTDRGWYCDNAGSVSFFPQADTCKVQSNRVFCDTMNSLTLPSEVSLCNTDIFNSKYD
                        CKIMTSKTDISSSVITSLGAIVSCYGKTKCTASNKNRGIIKTFSNGCDYVSNKGVDTV
                        SVGNTLYYVNKLEGKNLYVKGEPIINYYDPLVFPSDEFDASISQVNEKINQSLAFIRR
                        SDELLHNVNTGKSTTNIMITAIIIVIIVVLLSLIAIGLLLYCKAKNTPVTLSKDQLSG
                        INNIAFSK'''.replace('\n','').replace(' ','')

    ref = '''MELLIHRSSAIFLTLAINALYLTSSQNITEEFYQSTCSAVSRGY
                        LSALRTGWYTSVITIELSNIKETKCNGTDTKVKLIKQELDKYKNAVTELQLLMQNTPA
                        VNNRARREAPQYMNYTINTTKNLNVSISKKRKRRFLGFLLGVGSAIASGIAVSKVLHL
                        EGEVNKIKNALQLTNKAVVSLSNGVSVLTSRVLDLKNYINNQLLPMVNRQSCRISNIE
                        TVIEFQQKNSRLLEITREFSVNAGVTTPLSTYMLTNSELLSLINDMPITNDQKKLMSS
                        NVQIVRQQSYSIMSIIKEEVLAYVVQLPIYGVIDTPCWKLHTSPLCTTNIKEGSNICL
                        TRTDRGWYCDNAGSVSFFPQADTCKVQSNRVFCDTMNSLTLPSEVSLCNTDIFNSKYD
                        CKIMTSKTDISSSVITSLGAIVSCYGKTKCTASNKNRGIIKTFSNGCDYVSNKGVDTV
                        SVGNTLYYVNKLEGKNLYVKGEPIINYYDPLVFPSDEFDASISQVNEKINQSLAFIRR
                        SDELLHNVNTGKSTTNIMITAIIIVIIVVLLSLIAIGLLLYCKAKNTPVTLSKDQLSG
                        INNIAFSK'''.replace('\n','').replace(' ','')


sites = {"zero": list(range(62, 71)) + list(range(196, 211)),
            "one": [30, 32, 33, 39, 40, 42, 311, 312, 314, 343, 345, 348, 349, 376, 377, 379, 381, 382, 383, 390, 391, 467], #list(range(250, 259)) + list(range(261, 276)), 
            "two": list(range(254, 278)),
            "three": list(range(201, 214)) + list(range(429, 437)),
            "four": list(range(422, 438)),
            "five": list(range(131, 157)),
            "six": list(range(26, 37))
            }
#haoMuts = [12, 15, 19, 42, 103, 113, 116, 127, 173, 190, 191, 206, 209, 211, 220, 327, 389, 466, 467, 529]
antigenic_sites_data = {
    "Site Ø": {"ranges": [(62, 96), (195, 227)], "color": "#FFB6C1"},
    "Site I": {"ranges": [(27, 45), (312, 318), (378, 389)], "color": "#ADD8E6"},
    "Site II": {"ranges": [(254, 277)], "color": "#90EE90"},
    "Site III": {"ranges": [(46, 54), (301, 311), (345, 352), (367, 378)], "color": "#FFD700"},
    "Site IV": {"ranges": [(422, 471)], "color": "#DA70D6"},
    "Site V": {"ranges": [(55, 61), (146, 194), (287, 300)], "color": "#FAA460"},
    "P27 (FP)": {"ranges": [(110, 136)], "color": "#778899"},
#}
#antibody_sites_data = {
    "Nirsevimab": {"ranges": [(62, 69), (196, 212)], "color": "darkviolet", "linestyle":'-', "linewidth": 2.5, "label": "Nirsevimab", "label_fontsize": 12, "label_x_offset": 0.65},
    "Palivizumab": {"ranges": [(262, 275)], "color": "darkcyan", "linestyle":'--', "linewidth": 2.5, "label": "Palivizumab", "label_fontsize": 12, "label_x_offset": 0.65},
}

nextstrainDir = f"{baseDir}/NextstrainAnalysis"
combinedmeta = f'{nextstrainDir}/data/meta_{subtype}cc.csv'
meta = pd.read_csv(combinedmeta)

# New columns (to be updated by traversal)
meta['cluster'] = -1
meta['confidence'] = 0.0
meta['origin'] = ''
meta['clusterSize'] = 0
for site in sites:
    meta[f"as_{site}"] = 0 ## mutationscore
for locus in haoMuts:
    meta[f"asloc_{locus}"] = 0 ## mutationscore
meta['Region'] = [natDict.get(nat ,None) for nat in meta['Nationality']]

verbose = False
if subtype=='A':
    ali = f'{nextstrainDir}/data/alignmentA_3.fasta'
    aliF = f'{nextstrainDir}/data/alignmentA_3_Fprot.fasta'
    aliG = f'{nextstrainDir}/data/alignmentA_3_Gprot.fasta'
colors = pd.read_csv('../RSV/NextstrainAnalysis/config/colors.tsv', sep="\t", header=None)
colsfixed = [c.strip() for c in colors[2]]
colorDict = dict(zip(colors[1], colsfixed))
country_colors = {
    'philippines': '#1f77b4',  # Blue
    'france': '#ff7f0e',       # Orange
    'slovenia': '#2ca02c',     # Green
    'guangdong': '#d62728',    # Red
    'germany': '#9467bd',      # Purple
    'australia': '#8c564b',    # Brown
    'russia': '#e377c2',       # Pink
    'bangladesh': '#17becf'    # Cyan
}
colorDict.update(country_colors)
colDictUsed = {}
timestamps = []

# For iTOL export (piechart, collapsing)
itol_clades = []
collapsable_clades = []
collapsedCladeInfo = {}
allMutations = {}

fig, ax = plt.subplots(figsize=(12, 14))

qcountry = 'UAE'
def traverse(t, countryOrigInfo=None, targetClade=False, Fmuts=[]): 
    global subcladeCounter, timestamps
    global ax, verbose, colorDict, colDictUsed, itol_clades, collapsable_clades, collapsedCladeInfo
    global allMutations
    collapsable = True
    targetCladeRoot = False 
    subtreeInfo = []
    countryInfo = t['node_attrs']['country']
    Fmuts = list(Fmuts)  # Don't mutate parent list
    Fmuts += t['branch_attrs'].get('mutations', {}).get('F', [])
    country =  countryInfo['value']
    if not countryOrigInfo is None:
        countryOrig = countryOrigInfo['value']
    else: countryOrig = None
    if countryOrig != country and country == qcountry and countryInfo['confidence'][qcountry] > entryThreshold:
        putOrig = countryOrigInfo['value']
        confidence = countryOrigInfo['confidence'][putOrig]
        targetClade = True
        targetCladeRoot = True
    elif countryOrig == qcountry and country != qcountry:
        print('EXPORT:', country)

        
    # traverse children
    if 'children' in t: ## internal node
        childrenCollapsable = [True] * len(t['children'])
        colBranchStatsChildren = []
        for i, child in enumerate(t['children']):
            samples, childrenCollapsable[i], cbs = traverse(child, countryInfo, targetClade=targetClade, Fmuts=Fmuts)
            colBranchStatsChildren.append(cbs)
            if targetClade:
                subtreeInfo += samples
        collapsable = all(childrenCollapsable) # try to collapse all children, 
        #if t['name'] == 'NODE_0000257': pdb.set_trace()
        if not collapsable: # if just 1 not collapsable, collapse at least the ones that are
            #pdb.set_trace()
            for child, childCollapsable, cbs in zip(t['children'], childrenCollapsable, colBranchStatsChildren):
                if childCollapsable and 'children' in child:
                    #pdb.set_trace()
                    collapsedCladeInfo[child['name']] = cbs
                    collapsable_clades.append( child['name'])
        colBranchStats = Counter()
        for cbs in colBranchStatsChildren:
            colBranchStats += cbs

    else: ## leaf node
        nationality = t["node_attrs"].get('country', {'value':'?'})['value'] ## TODO: deal with metadata
        nodeName = f"{t['name']}"
        if country == qcountry:       
            subtreeInfo = [nodeName]
            collapsable = False ## UAE node, don't collapse the encompassing branches above
            
            mutations = collapseMutHist(Fmuts)
            mutsWRTnir = calcAntigenicDrift(mutations)
            #if len(mutations)<12:
            #    print(mutations,'->', mutsWRTnir)
            antigenStats = antigenicSiteMutScores(mutsWRTnir)
            #print(nodeName, mutations)
            allMutations[nodeName] = mutations
            #print(nodeName, antigenStats, len(mutations))
            for mutLoc in haoMuts:
                meta.loc[meta.strain==nodeName, f'asloc_{mutLoc}'] = int(mutLoc in mutations) # binary
            for site, mutScore in antigenStats.items():
                meta.loc[meta.strain==nodeName, f'as_{site}'] = mutScore # * 100
        elif targetClade:
            print(f'Exported to {nationality} ({nodeName})')
        colBranchStats = Counter([nationality])
        
    if targetCladeRoot and len(subtreeInfo) > 0:
        timestamp = t['node_attrs']['num_date']['value']
        timestamps.append(timestamp)
        coi = countryOrigInfo
        wedges =  np.array(list(coi['confidence'].values()))
        cols = [colorDict.get(country.lower(), '#A0A0A0') for country in coi['confidence'].keys()]
        td = dict(zip(coi['confidence'].keys(), cols))
        colDictUsed.update(td)
        Y = 0
        if len(subtreeInfo) < 17: Y = random.random()*15
        draw_pie(ax, wedges, timestamp, Y, colors=cols, size=max(len(subtreeInfo)/20, 0.1), alpha=max(confidence, 0.2))
        if True:
            #print(f'{putOrig} {100*confidence:.1f} {t["name"]}', end=' ')
            #print('-'*30, 'CLADE', subcladeCounter, '-'*50)
            #print(len(subtreeInfo), ','.join(subtreeInfo))
            #print('-'*100)
            ## adding new (calculated) columns to metadata 
            for node in subtreeInfo:
                if node in list(meta.strain):
                    meta.loc[meta.strain==node, 'cluster'] = subcladeCounter
                    meta.loc[meta.strain==node, 'origin'] = putOrig
                    meta.loc[meta.strain==node, 'confidence'] = confidence
                    meta.loc[meta.strain==node, 'clusterSize'] = len(subtreeInfo)
        # ----------- iTOL piechart export: collect data -------------
        itol_clades.append({
            'node': t['name'],
            'confidence': coi['confidence'].copy(),  # country: probability
            'clade_id': subcladeCounter
        })
        subcladeCounter += 1


    return subtreeInfo, collapsable, colBranchStats


# ---- File paths ----
rsv = {"A": '/mnt/Drive1/ahenschel/SambaShare/Flu_Madikay/RSV/NextstrainAnalysis/auspice/rsv_A_3b.json',
     "B": '/mnt/Drive1/ahenschel/SambaShare/Flu_Madikay/RSV/NextstrainAnalysis/auspice/rsv_B_5.json'}[subtype]

with open(rsv) as zf:
    tree = json.load(zf)

tree = tree['tree']
traverse(tree)

meta.to_csv(f'meta_{subtype}_antigenicStats.csv')
ax.hlines(y=0.5, xmin=min(timestamps) - 1, xmax=max(timestamps) + 1, color='black', linewidth=1, zorder=1)

ax.set_xlim(min(timestamps) - 1, max(timestamps) + 1)
legend_elements = [Patch(facecolor=col, label=country) for country, col in colDictUsed.items()]

ax.legend(handles=legend_elements, bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
ax.set_ylim(0, 1)
ax.set_yticks([])
ax.set_xticks(timestamps)
ax.set_xlabel("Timeline")
ax.set_title("Introductions to the UAE")
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig(f'intros2uae{subtype}.svg')
#plt.show()

# ---- Funtciont to write iTOL datasets
# pie chart dataset ----
def export_itol_piechart(itol_clades, colorDict, filename=f"itol_piechart_{subtype}.txt"):
    # Gather all countries
    all_countries = set()
    for entry in itol_clades:
        all_countries.update([c for c in entry['confidence'].keys()])
    all_countries = sorted(all_countries, key=lambda x: x.lower())
    field_colors = [colorDict.get(c.lower(), "#A0A0A0") for c in all_countries]
    # Build header
    header = [
        "DATASET_PIECHART",
        "SEPARATOR TAB",
        "DATASET_LABEL\tIntroductions to UAE",
        "COLOR\t#ff0000",
        "FIELD_LABELS\t" + "\t".join(all_countries),
        "FIELD_COLORS\t" + "\t".join(field_colors),
        "LEGEND_TITLE\tOrigin",
        "LEGEND_SHAPES\t" + "\t".join(["1"] * len(all_countries)),
        "LEGEND_COLORS\t" + "\t".join(field_colors),
        "LEGEND_LABELS\t" + "\t".join(all_countries),
        "",
        "DATA"
    ]
    # Data rows
    rows = []
    for entry in itol_clades:
        vals = [str(entry['confidence'].get(c, 0)) for c in all_countries]
        rows.append(f"{entry['node']}\t0\t20\t" + "\t".join(vals))
    with open(filename, "w") as out:
        out.write("\n".join(header + rows))
    print(f"iTOL piechart dataset written to {filename} ({len(itol_clades)} clades)")

def collapse_clades(intNodes2collapse, filename=f'collapse_{subtype}.txt'):
    header = '''COLLAPSE\nDATA'''
    with open(filename, 'w') as out:
        print(header, file=out)
        intNodes2collapse = [in2c.replace('_', '-') for in2c in intNodes2collapse]
        print('\n'.join(intNodes2collapse), file=out)
def prettyPrint(cbs):
    return '|'.join([f'{country}/{count}' for country, count in cbs.most_common(3)])
def writeCollapseCladeInfo(collapsedCladeInfo, filename=f'collapsedCladeLabels_{subtype}.txt'):
    header = '''LABELS\nSEPARATOR COMMA\nDATA'''
    with open(filename, 'w') as out:
        print(header, file=out)
        for node, cbs in collapsedCladeInfo.items():
            print(f'{node.replace("_","-")},{prettyPrint(cbs)}', file=out)

def branchColors(filename=f'branchcolors_{subtype}.txt'):
    header="""DATASET_RANGE
SEPARATOR COMMA
DATASET_LABEL,Clade Coloring
COLOR,#000000
RANGE_COVER,tree
DATA"""
    grays = ['#a0a0a0', '#c0c0c0']
    with open(filename, 'w') as out:
        print(header, file=out)
        for i, clade in enumerate(itol_clades):
            color = grays[i%(len(grays))]
            node = clade['node'].replace('_','-')
            print(f'{node},{node},{color},#ffffff,#000000,solid,1,UAE{i+1},#0000ff,1,italic', file=out)


def branchColors2(filename='nextcladeBranches.txt'):
    header="""DATASET_RANGE
SEPARATOR COMMA
DATASET_LABEL,Clade Coloring
COLOR,#000000
RANGE_COVER,tree
DATA"""
    with open(filename, 'w') as out:
        print(header, file=out)
        for i, clade in enumerate(itol_clades):
            color = '' ## to be fixed
            node = clade['node'].replace('_','-')
            print(f'{node},{node},{color},#ffffff,#000000,solid,1,UAE{i+1},#0000ff,1,italic', file=out)

def metaInfo(color_palette=color_palette, attr='clade'):
    filename = f'meta_{subtype}_{attr}.txt'
    uae = meta[meta.country=='UAE']
    header="""DATASET_COLORSTRIP
SEPARATOR COMMA
DATASET_LABEL,%s
COLOR,#ff00f0
STRIP_WIDTH,30
MARGIN,10
#BORDER_COLOR,#0000ff
SHOW_STRIP_LABELS,1
SIZE_FACTOR,14
LEGEND_TITLE,%s

""" % (attr,attr)    
    
    with open(filename, 'w') as out:
        print(header, file=out)
        print('LEGEND_SHAPES,' + ','.join(['1']*len(color_palette)), file=out)
        print('LEGEND_COLORS,' + ','.join(color_palette.values()), file=out)
        print('LEGEND_LABELS,' + ','.join(map(str, color_palette.keys())), file=out)
        print('DATA', file=out)
        for i, row in uae[['strain', attr]].dropna().iterrows():
            if row[attr] in color_palette:
                print(f"{row['strain']},{color_palette[row[attr]].lower()}", file=out)
    print('Completed:', filename)

def get_distinct_colors(items, cmap_name='tab20'): # Set 3
    n = len(items)
    cmap1 = cm.get_cmap(cmap_name, min(n, 20))  # Get colormap with n discrete colors    
    colors = [mcolors.to_hex(cmap1(i)) for i in range(min(20,n))]
    if n>20:
        cmap2 = cm.get_cmap('Paired', 10)  # Get colormap with n discrete colors
        colors += [mcolors.to_hex(cmap2(i)) for i in range(n-20)]
    return {item: color for item, color in zip(items, colors)}

export_itol_piechart(itol_clades, colorDict, filename=f"itol_piechart_{subtype}.txt")
collapse_clades(collapsable_clades)
writeCollapseCladeInfo(collapsedCladeInfo)
branchColors()

metaInfo(attr='clade')

uae = meta[meta.country=='UAE']
palettes = 'tab10 Pastel1 Set1 tab20b Set3 tab20 Set2 Set1'.split()

for p, attr in zip(palettes, 'Region Gender AgeGroup Facility ARI.SARI cluster clade city'.split()):
    items = sorted(uae[attr].dropna().unique())
    if attr=='AgeGroup':
        items = ['0-1', '1-5',  '5-10', '10-18', '18-35', '35-55', '55+']
        colors = '#FF0000,#FF7F00,#FFFF00,#00FF00,#0000FF,#8B00FF,#999999'.split(',')
        palette = dict(zip(items, colors))        
    elif attr=='city':
        colors = generate_blue_shades(n=len(items))
        palette = dict(zip(items, colors))         
    elif attr=='Region': # in shades of red
        colors = generate_blue_shades(hue=0, n=len(items))
        palette = dict(zip(items, colors)) 
    elif attr=='cluster':
        colors = gray_scale(len(items))
        palette = dict(zip(items, colors)) # * ((len(items)//10)+1)))     
    else:
        palette = get_distinct_colors(items, cmap_name=p)
    metaInfo(color_palette=palette, attr=attr)


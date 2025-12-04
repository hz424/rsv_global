# %%
# conda activate bio
from Bio import SeqIO
import numpy as np

def avg(l): 
    if not l: return 0
    return sum(l)/len(l)

seqDataDir = '/home/ahenschel/Flu_Madikay/RSV/tree/'
threshold = 5000 # rather inclusive
thresholdG = 20
reso = 5
nseqs = {'A': 99, 'B': 263} # Flu_Madikay/RSV/tree$ grep -c '^>' *aln_qc.fasta
genomeSize = 15254
genomeVis = genomeSize//reso
strain = 'B'

ref2fa = {'EPI_ISL_1653999': f'{seqDataDir}all_ref_EPI999.fasta', 
          'OK649754.1': f'{seqDataDir}all_ref_OK.fasta', ## the one with duplication
          'KT992094.1': f'{seqDataDir}All_RSV_{strain}_genomes_aln_qc.fasta'
          }
ref = 'EPI_ISL_1653999'
#ref = 'KT992094.1'
#ref = 'OK649754.1'

covDataDir = f'/mnt/Drive1/ahenschel/SambaShare/Flu_Madikay/HPC/Coverage_{ref}/'

totalCov = np.zeros((nseqs[strain], genomeVis))
totalCov1 = np.zeros((nseqs[strain], genomeSize))
covF = np.zeros((nseqs[strain], 1725))
covG = np.zeros((nseqs[strain], 897))

## input file (fasta files are assembled against different references)
fafile = ref2fa[ref]

## ouput files:
genomeFa = f'{seqDataDir}/all_RSV_{strain}_genomes_{ref}_masked.fasta'
ffa = f'{seqDataDir}/all_RSV_{strain}_genomes_{ref}_masked_F.fasta'
gfa = f'{seqDataDir}/all_RSV_{strain}_genomes_{ref}_masked_G.fasta'
# %%
filtered = 0
with open(genomeFa, 'w') as outfile, open(ffa, 'w') as outfileF, open(gfa, 'w') as outfileG:
    for row, rec in enumerate(SeqIO.parse(fafile, 'fasta')):
        if rec.id == "KT992094.1": continue
        if ref.startswith('EPI'): ## manually remove some unaligned seqs
            if rec.id[-5:] in ('12237', '02837', '12891', '11518'): continue
        if ref.startswith('OK') and len(rec.seq) != genomeSize: continue

        filename = f'{covDataDir}{rec.id}.s.f.bam.txt'
        mask = list('N'*genomeSize)
        covRow = []
        for line in open(filename):
            pos = int(line.split('\t')[1])-1
            cov = np.log(int(line.split('\t')[2])+1)
            mask[pos] = rec.seq[pos]
            totalCov1[row, pos] = cov
        seq = "".join(mask)
        totalCov[row] = [seq[i*reso:(i+1)*reso].count('N') for i in range(genomeVis)]
        #totalCov1[row] = covRow #[avg(covRow[i*reso:(i+1)*reso]) for i in range(genomeVis)]
        F = seq[5661:7386]
        G = seq[4688:5585]
        covF[row] = [int(nt=='N') for nt in F]
        covG[row] = [int(nt=='N') for nt in G]
        
        if seq.count('N') < threshold:
            #print(rec.id, seq.count('N'))
            print(f'>{rec.id}\n{seq}', file=outfile)
        else:
            filtered += 1
        if F.count('N') < thresholdG:
            print(f'>{rec.id}\n{seq}', file=outfileF)
        if G.count('N') < thresholdG:
            print(f'>{rec.id}\n{seq}', file=outfileG)
 
# %%
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, Normalize

masked_data = np.ma.masked_where(totalCov1 == 0, totalCov1)
cmap = cm.viridis  # or any colormap you like
cmap.set_bad(color='red')  

fig, axes = plt.subplots(5,1, figsize=(21,7))
x = range(0,genomeVis,200)
axes[0].set_xticks(x)  # Set ticks to match your x data
axes[0].set_xticklabels([f"{val * 10}" for val in x])
axes[0].matshow(totalCov)

axes[1].matshow(masked_data[:,:4000], cmap=cmap)
axes[2].matshow(masked_data[:,4000:8000], cmap=cmap)
axes[3].matshow(masked_data[:,8000:12000], cmap=cmap)
axes[4].matshow(masked_data[:,12000:], cmap=cmap)

#axes[2].matshow(covF)
#axes[3].matshow(covG)
plt.tight_layout()
plt.savefig(f'coverage_rsv{strain}_{ref}.svg')

# %%


plt.figure(figsize=(14,9))
plt.matshow(totalCov)
#axes[1].matshow(covF)
#axes[2].matshow(covG)
plt.tight_layout()
# %%

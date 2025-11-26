import pandas as pd
from scipy import stats
import numpy as np

# Leggi il file tabulato
megahit = pd.read_csv("megahit_contigs_distribution.txt", sep="\t", names=["contig", "length"]).loc[:,"length"]
spades = pd.read_csv("spades_contigs_distribution.txt", sep="\t", names=["contig", "length"]).loc[:,"length"]
trinity = pd.read_csv("trinity_contigs_distribution.txt", sep="\t", names=["contig", "length"]).loc[:,"length"]

megahit = pd.to_numeric(megahit)
spades = pd.to_numeric(spades)
trinity = pd.to_numeric(trinity)

mega_spa = stats.mannwhitneyu(megahit, spades, alternative="greater").pvalue
mega_tri = stats.mannwhitneyu(megahit, trinity, alternative="greater").pvalue
spa_tri = stats.mannwhitneyu(spades, trinity, alternative="greater").pvalue

with open("stat_testing.log", "w") as f:
	f.write(f"""Using Mann-Whitney U test to test if assemblers create contigs of different length.
	Megahit:
		Contigs assembled: {megahit.count()}
		Average length:{megahit.mean()}
		Standard deviation:{megahit.std()}
	Spades:
		Contigs assembled: {spades.count()}
		Average length: {spades.mean()}
		Standard deviation: {spades.std()}
	Trinity:
		Contigs assembled: {trinity.count()}
		Average length: {trinity.mean()}
		Standard deviation: {trinity.std()}\n
	""")
	f.write(f"\n")
	f.write(f"""Result for wilcoxon ranked test between megahit and spades: 
			Megahit assembles, on average, longer contigs than spades with p-value = {mega_spa:.15e}\n""")
	f.write(f"""Result for wilcoxon ranked test between megahit and trinity:
			Megahit assembles, on average, longer contigs than trinity with p-value = {mega_tri:.15e}\n""")
	f.write(f"""Result for wilcoxon ranked test between spades and trinity: 
			Spades assembles, on average, longer contigs than trinity with p-value = {spa_tri:.15e}\n""")

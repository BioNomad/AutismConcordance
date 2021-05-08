# AutismConcordance

> Here we provide the repository for the data analysis performed by in our work "Personalized Perturbation Profiles Reveal Concordance between Autism Blood Transcriptome Datasets" <sup>1</sup>. In this analysis we take publically available blood transcriptomic datasets and attempt to establish concordance between them. Historically this has been a difficult task given the heterogeneity of the disease. We leverage personalized pertubation profiles to establish concordance <sup>1</sup>:
<p align="center">
    <img src="https://user-images.githubusercontent.com/59677194/117555115-6edca100-b02a-11eb-8d3a-bc8c413863e4.png" width="500" height="650" />
</p>
> Here a gene is considered perturbed in a study if its expression value z-score is greater than 2.5 or less than -2.5. This gene then gets added to each study's "list" of perturbed genes. We then compare these lists and found a common set of perturbed genes. Moreover, we found that this common pool of perturbed genes contains the following poorly characterized genes: C18orf25, C15orf39, C1orf109, C1orf43, C19orf12, C6orf106, C3orf58, C19orf53, C17orf80, C4orf33, C21orf2, C10orf2, C1orf162, C10orf25 and C10orf90 <sup>1</sup>. We then use differential correlation analysis to identify how these genes might be relevant to autism molecular mechanisms <sup>1</sup>:
<p align="center">
    <img src="https://user-images.githubusercontent.com/59677194/117555112-64baa280-b02a-11eb-9a23-ed1c74af4baa.png" width="500" height="650" />   
</p>
To confirm these correlations we identify those seen in TCGA Low Grade Glioma data. These correlations were explored with the text mining tool Chilibot to reveal interesting connections to DNA damage, ubiquitination, R-loops, autophagy, and mitochondrial damage <sup>1</sup>. 


# References

> 1. Laird, J. & Maertens, A. Personalized Perturbation Profiles Reveal Concordance between Autism Blood Transcriptome Datasets. bioRxiv 2021.01.25.427953 (2021) doi:10.1101/2021.01.25.427953.

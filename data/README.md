
<head>
<style>
h1 {font-size: 16pt;
    font-family: arial, helvetica, sans-serif;
    }
h3 {font-size: 14 pt;
    font-family: arial, helvetica, sans-serif;
}
body { background-color: white ; font-size: 11pt;
    font-family: arial, helvetica, sans-serif;
 }
td {font-size: 11pt; text-align:left;
    font-family: arial, helvetica, sans-serif;
}
th {font-size: 11pt;
    font-family: arial, helvetica, sans-serif;
}
</style>
</head>

## FASTA scoring matrices in `data/`

The `.mat` files in this directory show the similarity scores
associated with the various `-s abbrev` scoring options for the
`fasta36` programs.  Thus, `-s BP62` uses the *BLOSUM62* matrix, which
is provided in `blosum62.mat`.  There are four families of matrices provided in this directory:

1. The `pam120.mat` and `pam250.mat` matrices derived from the original Dayhoff (1978) PAM data (Dayhoff et al. (1978), Atlas of Protein Sequence and Structure, Vol. 5 suppl 3 pp. 345-352; National Biomedical Research Foundation, Silver Spring, MD

2. The BLOSUM series of matrices (`blosum50.mat`, `blosum62.mat`, etc., Henikoff and Henikoff (1992), PNAS 89:10915-10919)

3. The MD series of matrices (`md_10.mat`, `md_20.mat`, `md_40.mat`) described by Jones, Taylor, and Thornton (1992) Comp. Appl. Biosci 8:275-282.

4. A series of matrices developed using the VTML model (Mueller and Vingron (2002) Mol. Biol. Evol. 19:8-13): `VTML_10.mat`, `VTML_20.mat`, etc.

5. In addition, the Optima5 matrix (Kahn and Goldstein (2002) Proteins 48:367-376) is provided (`-s OPT5`).

---

The mapping of `-s abbrev` to `data/*.mat` file is found in
`src/upam.h` file in the definition of the `struct std_pam_str std_pams[]`
file:


| Abbrev. |  .mat file  |
| ------- | ------------ |
| `VT10`  | `vtml_10.mat`|
| `P10`   | `md_10.mat`  |
| `M10`   | `md_10.mat`  |
| `MD10`  | `md_10.mat`  |
| `VT20`  | `vtml_20.mat`|
| `P20`   | `md_20.mat`  |
| `M20`   | `md_20.mat`  |
| `MD20`  | `md_20.mat`  |
| `VT40`  | `vtml_40.mat`|
| `P40`   | `md_40.mat`  |
| `MD40`  | `md_40.mat`  |
| `VT80`  | `vtml_80.mat`|
| `BL80`  | `blosum80.mat`  |
| `VT120` | `vtml_120.mat` |
| `P120`  | `md_10.mat`  |
| `BL62`  | `blosum62.mat` with -8/-1 gap penalties |
| `BP62`  | `blosum62.mat` with -11/-1 gap penalties |
| `VT160` | `vtml_160.mat` |
| `BL50`  | `blosum50.mat` |
| `OPT5`  | `optima5.mat` |
| `PAM250` | `pam250.mat` |


---

In addition, several other scoring matrices are provided, including:

| `.mat` |     |
| ----   | ----|
| `idn_aa.mat` | identity matrix for protein searches |
| `dna.mat`    | +5/-4 matrix for DNA searches |
| `rna.mat`    | +5/-4 matrix for RNA searches with +2 for G:A/TU:C |


# iSFAC_code
Code used in the preparation of the Manuscript on iSFAC modelling
- Rscat contains a python script to compute Rscat (Yonekura et al, IUCrJ (2018),
  5, 348-353 https://doi.org/10.1107/S2052252518005237)

- pearson_coefficient contains a python script to compute the Pearson
  correlation coefficient between sets of partial charges

- tutorials contains a video on how to create and use SFAC command for SHELXL,
  both for ionic and neutral atoms. Unfortunately, the volume of this video is
  very low

- cubemaps contains code for reading files and computing their Pearson
  correlation coefficient. Cube maps can be on different grids. Atoms will be
  superimposed with the Kabsch algorithm. This requires that the N atoms of the
  reference map (excluding hydrogen atoms) are identical and in the same order
  as the first N atoms of the reference map.

- fcf2cube converts the SHELXL FCF-file to a CUBE file, covering the atoms
  within the corresponding RES file. It should avoid 'dark patches' on the map
  shown by the program VMD, that are observed when coordinates lie outside the
  crystallographic unit cell.


The following code follows from Nasirudin et al. (2019) where we use a high resolution simulation (in our case one ran using cubeP3M)
in order to predict the collapse fractions in a lower resolution simulation.

If the hgh-resolution simulation was ran using cubeP3M you can use cond_mass_funct_david.f90 which will generate files containing the 
collapse farctions and density in each coarened grid. The program can also generate the collpase fraction in different mass bins.

Once you have the binned data, you will also need the density field of the low res simulation. Once you have these two files you can run
fcoll_implementation.py which will generate files of the scatter implemented collapse fractions.


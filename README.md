# I. Protocol for Running Pareto Scans and Plotting Results

### Establish remote connection to corin from your computer
- If using Windows, open a new PuTTY session, enter shivi@149.151.162.35 in the hostname field and click 'Open'.  If prompted with a question, select "Yes".  Enter password when prompted.
- If using a Mac, you can follow this video tutorial: https://www.youtube.com/watch?v=DpgQe_j371E&feature=youtu.be 

### Clone github repository (first time only)
- To clone github repository, type

`git clone https://github.com/jayfoleyiv/LightBulb`

- To change into the LightBulb directory, type
`cd LightBulb`

## The following instructions use the noble metal-oxide system as an example (see system 3 under II. Notes on Different Systems)

- Change directories to the MetalOxide folder (from within the LightBulb folder)

`cd MetalOxide`

- To view the contents of this folder, use the 'ls' command (ls for list)

`ls`

- At a minimum, you should see the following files/folders after typing ls:

`DIEL  input.txt  plot_pareto.gnu  plot_spectra.gnu  Scan.c  Scan.exe  SPECTRA`

- Scan.exe is a program that will scan through a number of structural and material parameters, identify the set of Pareto optimal structures, and compute the spectral efficiency, useful power density, and thermal emission spectrum for all Pareto optimal structures

- input.txt is an input file that specifies the material and structural parameters that Scan.exe will use.  Before running Scan.exe, open input.txt and verify it contains parameters relevant to the system you want to study.  You can use the text editor 'vim' to open input.txt:

`vim input.txt`

- <b> NOTE: </b> vim is a text editor, for a tutorial on basic usage of vim following [this link.](http://www.openvim.com/)

- Once input.txt is open, you will see a variety of parameter labels (words) and parameter values (numbers).  An example of input.txt is provided in Section <b> II. Input File </b> for reference, and the explanation of the labels and values are commented in this readme for clarity.  Comments are in parenthesis and should <b> NOT </b> appear in the actual file input.txt  
Let's say you've decided to perform a scan of structures where the metal is Rhodium (Rh), and have changed the input file appropriately such that the file prefix is <i> Rh_test </i> (see <b> II. Input File </b> for more information about how to change the input file)

- Run the program Scan.exe and direct the output to a file called 'Rh_Scan.txt'

`./Scan.exe input.txt >& Rh_Scan.txt &`

<b> NOTE: </b> The program will scan millions of structures, so it will take hours to finish.  If you start a calculation in the morning, check back in the evening.  If you start it in the evening, check back in the morning.

- To plot the Pareto front which will be stored in a file called <i> Rh_test_Pareto.txt </i>, first open up the file <i> plot_pareto.gnu </i>

`vim plot_pareto.gnu`

- Edit the contents of <i> plot_pareto.gnu </i> to match the following:

```gnuplot 
#!/usr/bin/gnuplot

set terminal postscript enhanced color 'Helvetica' 25
set output 'Rh_test_Pareto.eps'
set xlabel 'Spectral Efficiency'
set ylabel 'Sectral Density (W/m^2)'
set pointsize 3
plot 'Rh_test_Pareto.txt' u 6:7 w p
```

- Execute the gnuplot script you just edited by typing

`./plot_pareto.gnu`

- There will typically be many points on the Pareto front, to plot the thermal emission spectra of the structure corresponding to the 5th point on the Pareto front, edit the content of <i> plot_spectra.gnu </i> to match the following:

```gnuplot 
#!/usr/bin/gnuplot

set terminal postscript enhanced color 'Helvetica' 25
set output 'Rh_test5_spectra.eps'
set xlabel 'Wavelength (nm)'
set ylabel 'Thermal Emission (W/ m^2 / nm / sr )'
set pointsize 3
plot 'Rh_test5.txt' u 1:4 w l lw 4 title 'TE Spectrum', \
'Rh_test5.txt' u 1:5 w l lw 4 title 'Blackbody Spectrum'
```

- Execute the gnuplot script you just edited by typing

`./plot_spectra.gnu`

- Once you have plotted all quantities you wish to analyze, push your changes to github:

`git init`

`git add .`

`git commit -m 'Added my project'`

`git push`

- Enter Username when prompted!

## II. Input File
- Note: Comments about the meaning of the input parameters are in parenthesis and should <b> NOT </b> appear in the actual file input.txt

>Nmin_Nmax   (minimum and maximum number of layers in the structure) 
>
>5 28        (scan will consider structures with number of layers between 5 and 28)

>d1min_d1max (minimum and maximum thickness of low refractive index layer) 
>
>0.1  0.3    (scan will consider low refractive index layers between 0.1 and 0.3 microns thick)
 
>d2min_d2max  (minimum and maximum thickness of high refractive index layer) 
>
>0.1  0.3     (scan will consider high refractive index layers between 0.1 and 0.3 microns thick)

>vfmin_vfmax  (minimum and maximum volume fraction of alloy layer)
>
>0.0 1.0      (scan will consider alloys with volume fractions between 0 and 100% of the metal in the oxide)

>Tmin_Tmax    (minimum and maximum temperature of the structures)
>
>1000 1700    (Scan will consider performance of structures when they are between 1000 and 1700 K)

>FilePrefix   (Controls what the name of various output files will be)
>
>Rh_Alumina   (All output files created by Scan.exe will begin with the word 'Rh_Alumina')

>AbsorberFileName (Specifies where Scan.exe should read the data that defines the metal material)
>
>DIEL/Rh_Spline.txt (Scan.exe will read the metal material data from the file Rh_Spline.txt in the DIEL folder... this is data for the metal Rhodium (Rh))


# III. Notes on Different Systems
(1) oxide-oxide layers (ex. Fe3O4-Al2O3, or Fe3O4-Fe2O3) on top of BR

    Fe3O4 is absorbing and has a high melting point.
    Oxidation states has to be well controlled (Fe2O3: transparent).
    There is an intermediate growth condition where both Fe2O3 and Fe3O4 formed.

(2) nitride-nitride layer buried under BR (ex. HfN-AlN)

    AlN, TiN, and HfN are resistive to oxidation when buried under thick oxide layers (hundreds of nanometer). 
    Nitrides may be a better barrier layer for W diffusion outward to oxide BR layers from W substrate.
    However, the effect of coupling to oxide BR is less dramatic when the weak absorber is below BR compared to when it’s above the BR.
    Precise control over stoichiometry at high temperature is needed. (Metal rich: absorbing, N rich: transparent)

(3) noble metal-oxide layers (ex. Ru-Al2O3) on top of BR

    Ru has a high melting point and similar optical properties to W
    RuO2 also has a higher melting point than other noble metal oxides (not sufficient though).

Attached slides contain my preliminary literature surveys on thermal stability and optical properties of the materials mentioned above. FYI, I also attached the refractive indices data obtained from the website: https://refractiveindex.info. So, my main question to you is which materials systems you expect to work best in terms of optical properties. For oxide BR, we agreed to use Al2O3 for low n, and ZrO2 or HfO2 for high n for better thermal stability. I'll measure the refractive indices of these oxides grown by PEALD for next runs of simulations.

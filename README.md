This is the coalescence+fragmentation hadronization code. When you use this code, please cite our papers:
       W.~Zhao, etc. al.[arXiv:2103.14657 [hep-ph]].
       W.~Zhao, etc. al.Phys. Rev. Lett.125, (2020) no.7, 072301.
       
       
######### Coalescence code: #############
The format for the input thermal parton for coalescence is:
     number_of_thermal_parton
     id, px, py, pz, energy, x, y, z, t

The format for the input shower parton for coalescence is:
     number_of_shower_parton
     id, px, py, pz, mass, x, y, z, t

The folder, Coal, contains the coalescence code.
You compile the coalescence code by:
        make
Then run the coalescence code by:
        ./main THE/PATH/TO/INPUT/PARTON
It will output the coaleced_hadron.dat that contains the coaleced hadrons. The file remnant_jet_parton.dat is the remnant partons that will hadronized by fragmentation.

######## Install Pythia8  ###########
     - Create a folder for installing Pythia8
     - Got to the created folder
     - wget https://pythia.org/download/pythia83/pythia8306.tar
     - tar -xvzf pythia8306.tgz
     - Go to the extracted folder 
          - cd pythia8306
     - Run make command
          - ./configure
          - make
          
######### Fragmentatoin code: #############
You first run the code in the folder, to_frag2_with_thth.
Compile the code:
        g++ -o tofrag main-to-urqmd-all.cc
Move the remnant_jet_parton.dat from Coal to the to_frag2_with_thth. Then run:
        ./tofrag
It will return the file jet_parton1.dat that for the string fragmentation.

cp -r main2000.cc TO/PYTHIA/EXAMPLE/FILE
cd TO/PYTHIA/EXAMPLE/FILE
First compile the pythia8 fragmentation code, main2000.cc by:
        make main2000
Then run the code:
        ./main2000 number_of_event THE/PATH/TO/hadrons_frag1.dat THE/PATH/TO/jet_parton1.dat 
The file hadrons_frag1.dat the hadrons from the string fragmentation.


######### Combine the coaleced hadron and fragmented hadron code: #############
The folder, to_UrQMD, is the code that combines the coaleced hadrons and fragmented hadron for the UrQMD hadronic evolution.
Compile the combine code:
        g++ -o tourqmd main-to-urqmd-all.cc
Then move the file coaleced_hadron.dat in the Coal into the folder to_UrQMD. Move the file hadrons_frag1.dat from the fragmentation into to_UrQMD folder. Then run:
        ./tourqmd number_of_event 
It will output the file oscar.dat that input for the UrQMD.

The format for the oscar.dat is:
     event_id    number_of_hadron
     id, px, py, pz, energy, mass, x, y, z, t

# Michael Joseph S. Lee's Master's Thesis Repository

## Description
This is a repository for Michael Joseph S. Lee's master's thesis.

***Note that this repository by itself does not contain everything—some code is contained in other repositories.  See the "Contents" section of this README for more details.***

## Contents
- **``MJSL Thesis Final Draft 7-23-2023.pdf``** - A PDF of my master's thesis
- **``Code/``** - Contains the code I used in my thesis
  - **``Code/DORT Background/``** - Contains the MATLAB scripts that I used to make the figures demonstrating DORT in the "DORT Background" chapter of my thesis
  - **``Code/MATLAB Simulations/``** - Contains the MATLAB scripts that I used to do the experiments in the "MATLAB Simulations" chapter of my thesis
  - **``Code/Imaging Through Aberrating Media Simulation/``** - Contains a link (in ``Code/Imaging Through Aberrating Media Simulation/README.md``) to the ``SW-DORT-Research-Code`` repository, which contains the code that I used to perform the experiment in the "Imaging Through Aberrating Media Simulation" chapter of my thesis
    - The ``SW-DORT-Research-Code`` repository can be found at https://github.com/mich-lee/SW-DORT-Research-Code/ .  Make sure to use Commit 463d8a3debcebdee6b7c5fddc2cfaea3828dac10 of that repository as that corresponds to what was used in this thesis.
    - The ``SW-DORT-Research-Code`` repository requires the repository at https://github.com/mich-lee/holofork
      - Download ``holofork`` repository (specifically Commit a27e3fe6990e81a4642bae3e6511e9d6385620e1 of https://github.com/mich-lee/holofork) and have it as the ``holotorch-lib`` folder in the downloaded ``SW-DORT-Research-Code`` repository
    - The ``THESIS EXPERIMENTS`` folder in the ``SW-DORT-Research-Code`` repository contains code that I used to generate and analyze data for the experiment in the "Imaging Through Aberrating Media Simulation" chapter of my thesis
      - ``SW-DORT-Research-Code/THESIS EXPERIMENTS/Generate_Data (aberrating layer 7mm in front of scatterer plane).py`` - Generates data
        - This script will generate a file containing data.
      - ``SW-DORT-Research-Code/THESIS EXPERIMENTS/Analyze_Data (THESIS).py`` - Analyzes the data generated by ``SW-DORT-Research-Code/THESIS EXPERIMENTS/Generate_Data (aberrating layer 7mm in front of scatterer plane).py``
        - This script contains code for plotting.  The plots generated by this script were used in the "Imaging Through Aberrating Media Simulation" chapter of my thesis.
        - When using this script, set the "dataFilePath" variable in this script to point to the data generated by the ``SW-DORT-Research-Code/THESIS EXPERIMENTS/Generate_Data (aberrating layer 7mm in front of scatterer plane).py`` script
- **``Conference Presentation/``** - Contains the conference presentation for a COSI 2023 paper that was accepted
  - **``Conference Presentation/COSI 2023 Presentation (Smaller Size).mp4``** - A video of my COSI 2023 paper presentation
  - **``Conference Presentation/COSI 2023 Presentation.pptx``** - The PowerPoint presentation for my COSI 2023 paper presentation
    - The presenters notes contain a transcript of what was said during the COSI 2023 paper presentation.
    - This presentation was derived from my thesis defense presentation at **``Thesis Defense Presentation/Thesis Defense Presentation.pptx``**.  Media files for my **``Conference Presentation/COSI 2023 Presentation.pptx``** presentation can be found in **``Thesis Defense Presentation/Media/``**.
- **``Thesis Defense Presentation/``** - Contains my thesis defense presentation
  - **``Thesis Defense Presentation/Thesis Defense Presentation.pptx``** - The PowerPoint presentation for my thesis defense presentation
    - The presenters notes contain a transcript of what was said during my thesis defense presentation.
    - Media files for my **``Thesis Defense Presentation/Thesis Defense Presentation.pptx``** presentation can be found in **``Thesis Defense Presentation/Media/``**.
  - **``Thesis Defense Presentation/Media/``** - Contains media files for my **``Thesis Defense Presentation/Thesis Defense Presentation.pptx``** presentation

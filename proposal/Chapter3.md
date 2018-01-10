---
output: pdf_document
fontsize: 11pt
#bibliography: ch3.bib #C:/Users/Julie/Box Sync/R/Proposal/ch3.bib
#csl: the-auk.csl #C:/Users/Julie/Box Sync/R/Proposal/the-auk.csl
---


#Phylogenetic versus morphological signal in movement characteristics of avian taxa 
##Introduction
Closely related species typically share many morphological, ecological, and behavioral characteristics, which are inherited from common ancestors. Shared traits among related species [i.e. phylogenetic signals, @blomberg2002tempo] have been found in a wide range of taxa and for a diversity of inherited traits. Phylogenetic signals have been found for sensitivity to pesticides in amphibians [@chiari2015phylogenetic], climatic niches of hylid frogs in the neotropics [@wiens2006evolutionary], growth and habitat use in rainforest plants [@chazdon2003tropical], and shifts in flowering time in plant communities of the Northern Hemisphere [@rafferty2017global]. 

Not all traits are conserved among related species. For example, characteristics related to the structure of vocalizations showed phylogenetic signals in herons, but frequency characteristics did not show phylogenetic signals and instead was predicted by the herons’ habitat [@mccracken1997avian]. Likewise, @kamilar2013phylogenetic showed phylogenetic signals of some morphological and behavioral characteristics of primates but no support for other behavioral and ecological characteristics. Traits that are under selection by the environment are likely to show low phylogenetic signals, including: climatic and habitat niches of European birds [@pearman2014phylogenetic], thermal tolerance of beetles [@garcia2016limited], and flowering response to snowmelt and temperature of a subalpine plant community [@caradonna2015phenological]. A major question in ecology, therefore, is which traits are due to shared inheritance from a common ancestor? And which traits are in response to an organism’s environment and ecology?

One trait that has not been well tested for phylogenetic signal is movement. Movement is determined, in part, by an organism’s movement ability (I.e. wing morphology and flight behavior), internal state (i.e. age/experience, migratory behavior, breeding status), and its environment (i.e. interspecific interactions, resources, habitat, uplift, weather; Nathan et al. 2008). It is predicted that shared, inherited traits that affect movement (e.g. morphology, flight behavior, migratory behavior) will show a phylogenetic signal where closely related species will move more similarly than distantly related species. In contrast, I predict that traits related to environment (e.g. resources, habitat, uplift, weather) or traits that change throughout an individual’s life (e.g. age, breeding status, interspecific interactions) will show an ecological signal. 
To determine which traits are due to a shared inheritance, I will calculate the mean value for several movement characteristics from track data of individual birds and morphological characteristics from specimens and field data. The variances of each trait will be compared to the relatedness of individuals from a phylogenetic tree. When trait variance is low compared to relatedness, there is a phylogenetic signal of that trait. 

Species may be more similar across multiple characteristics, instead of a single trait. Furthermore, similarity may be conserved across different levels of taxonomic organization. Therefore, it is necessary to determine at what level is movement predicted by taxonomy. To test this question, I will create a tree based on shared movement characteristics and determine at what level do related individuals cluster together. The null hypothesis, movement is individualistic, would be supported if individuals of a species did not cluster together. Alternatively, movement could be species specific, which would be supported if individuals cluster together at the species level. Lastly, movement could be shared at a higher level of taxonomy (i.e. genus or family); this would be supported if individuals of a genus or family cluster together.

##Methods
###Study Species
Due to the cost and the weight of telemetry units being proportional to the battery efficiency and sampling frequency, our dataset is biased towards species that are larger and can support the weight of a telemetry unit (3% rule). The 36 study species are listed in Table 1.

###Subsampling Data
To isolate similarities in movement that are due to 1) physical ability, 2) behavior, and 3) environment, I will remove location data that are associated with migration and nesting. To standardize movement data of birds with differing migratory behaviors, I will remove migratory tracks that were determined using a fourier analysis. Because not all birds were active breeders, I will remove nesting periods from tracks to remove the central foraging signal from the data. In short, I will only use 'winter' or non-breeding foraging movements of individual birds.

To standardize movement data that is collected at different frequencies, I will only use tracks that were sampled at hourly or finer temporal scales (I.e. >240 points in a 10 day interval). All tracks will be resampled to hourly intervals (xts R package). To standardize different sampling efforts (i.e. number of days of tracking), I will subsampled tracks of each individual to one continuous 10 day segment. Species that have few useable tracks (i.e. 3) will be removed from analyses. 

###Movement Characteristics Data
I will then determine each individual’s movement characteristic for this 10 day track. I will calculate mean and standard deviation of hourly displacement, hourly relative angle, and daily activity. ‘Daily activity’ is defined as the number of hours displacement was greater than the margin of error of the telemetry unit (i.e. 18 m).  I will also calculate, for the 10 day track, 1) the 90% minimum convex polygon, 2) the maximum displacement from the first location, and 3) the number of hours for first passage time of a 500 m radius (individuals that did not leave a 500 m radius will be assigned the maximum number of hours, i.e. 240).

###Phylogenetic and Morphological Data
I will collect data on wing morphology from one of the following sources: 1) author field data, 2) literature, 3) museum specimens. Wing characteristics I will test for phylogenetic signals are: body mass, wing loading, aspect ratio, wing length, and wing chord.
Avian phylogeny will be reconstructed from the latest avian phylogeny, created using next-generation sequencing [@prum2015comprehensive].

###Analyses
I will test for phylogenetic signals of mean value for each movement characteristic at the species level. This is because some methods for detecting phylogenetic signals cannot handle polytomies.

To test for the similarity of movement by taxonomic level, I will create a dissimilarity matrix (R package vegan, @vegan) for all of the individuals to be used in hierarchical clustering. To control for uneven sampling of species, which increases variance and decreases species level clustering, I will calculate a threshold of dissimilarity by number of individuals based on the expected variance for each sample size. 	

##Expected Results and Discussion
###Hierarchical clustering
I have 37 species and a total of 680 individuals to use in my analysis. These species represent several avian guilds. At our sampling resolution (hourly over 10 days), I generally expect movement characteristics to cluster due to ecology and guild. For example, waterfowl would cluster together (*Anas* with *Anser*, etc) because of the similarity in their habitats, flight mode, and body size. Finer scale differences that are related to their ecological niches should not be captured in our sampling. 

Likewise, I expect several vulture species to cluster together because they are large soaring species that search on the wing for carrion. However, species in my dataset are from both the old and new world vultures, so vultures may cluster either by ecological similarities (e.g. Egyptian vulture with new world vultures and other old world as a separate cluster) or by relatedness (i.e. new world species separate from old world species).

I also expect to find low clustering within many species due to environmental covariates affecting individual movement. I tested this using a hierarchical clustering of movement characteristics on a subset of the species in my dataset. I found individuals of two species clustered as at the species level: Oilbirds and Knob-billed ducks (Figure 1). No other species had all individuals clustered together. These results suggest movement characteristics may be more stereotypical in Oilbirds and Knob-billed ducks than in the other six species tested, which may be more strongly affected by the environment. 

###Phylogenetic signal
I expect to find low phylogenetic signals for many movement characteristics, due to some species and individuals being more strongly influenced by their environment or internal state than others. I expect wing morphology to show a stronger phylogenetic signal, but the strength of this relationship will be low due to convergence.

Wing loading and aspect ratio are the two primary measures of wing morphology. Together, they can inform about the ecology of a species. Some studies have looked at greater detail than these measures, to describe how small changes in wing morphology correspond to ecological and behavioral differences among related species, such as pelagic birds [@brewer2007wing].

I expect that wing morphology will explain maximum flight speed or distance for some species - those that require frequent fast or long distance flights. It is known that wings are under strong selective pressures. Soaring birds, for example, require long and wide wings with low wing loadings to soar efficiently [@hedenstrom1993migration]. Soaring birds are able to cover large distances at low energetic costs and often need to cover large distances to look for food or because of the network of updraft availability. Soaring birds, therefore, should have longer flight distances than other birds, and possibly greater flight speeds. 

For others species that deviate from the phylogenetic signal, I will need to follow up this analysis and determine if there are any life history, behavior, or environmental traits that would require birds to have a wing morphology that is different from their flight and movement patterns. 

##Justification
To my knowledge, no studies have looked for phylogenetic signals of movement, specifically of avian movement. First, it is interesting to learn if any movement characteristics are conserved, phylogenetically. If so, how conserved are they? Specifically, to what level of taxonomic organization? 

There are clear applications of this approach. First, many movement studies use data from several individuals to estimate the movement and ecological characteristics of a population or a species. It would be beneficial to know whether the study of a closely related species could be applied to an understudied species. Especially because logistical issues restrict data collection and sampling, including: cost, time, and success and precision of devices. Additionally, the biology can further restrict potential data collection: weight and battery life of units, difficulty and sensitivity of organisms to trapping, mortality of individuals, etc. Therefore, if any characteristics are phylogenetically conserved, it may be possible to estimate movement characteristics from a related species.


\newpage
##Tables
Table 1. Species (n=37) and number of individuals (n=680) available to be used in our analyses. Number of studies indicates whether data for species are available from one or more unique datasets.

|Species |Common Name | Number of studies | Number of individuals|
|------|-------|---|---|
|Anas crecca |Eurasian Teal| 1 | 4 |
|Anas penelope|Eurasian wigeon|1|1|
|Anas platyrhynchos|Mallard|2|16|
|Anas strepera|Gadwall|1|6|
|Anser cygnoides|Swan goose|1|15|
|Anser indicus|Bar-headed goose|2|5|
|Aptenodytes patagonicus|King penguin|1|9|
|Ardea alba|Great Egret|1|8
|Ardea herodias|Great Blue Heron|1|3|
|Calonectris diomedea|Cory's shearwater|1|15|
|Cathartes aura|Turkey Vulture|2|36|
|Ciconia ciconia|White Stork|5|70|
|Columba livia|Homing Pigeons|3|72|
|Coragyps atratus|Black Vulture|1|2|
|Creagrus furcatus|Swallow-tailed gull|1|41|
|Fregata magnificens|Magnificent Frigate Bird|1|9|
|Grus grus|Common Crane|1|20|
|Grus nigricollis|Black necked crane|1|4|
|Gyps africanus|White Backed Vulture|2|14|
|Gyps coprotheres|Cape Vulture|1|4|
|Gyps fulvus|Griffon Vulture|1|6|
|Haliaeetus leucocephalus|Bald Eagle|2|5|
|Morus capensis|Cape Gannet|2|20|
|Neophron percnopterus|Egyptian Vulture|1|2|
|Pandion haliaetus|Osprey|2|17|
|Phaethon aethereus|Red-billed tropicbird|1|19|
|Pteroglossus torquatus|Collared Aracari|1|3|
|Pygoscelis adeliae|Yan Penguin|1|6|
|Ramphastos sulfuratus|Keel-billed toucan|1|4|
|Ramphastos swainsonii|Chestnut-mandibled toucan|1|2|
|Rissa tridactyla|Kittiwake|1|10|
|Sarkidiornis melanotos|Knob-billed duck|1|7
|Steatornis caripensis|Oilbird|1|6|
|Sula dactylatra|Masked Booby|3|145|
|Sula leucogaster|Brown Booby|1|8|
|Sula sula|Red-footed Booby|1|62|
|Torgos tracheliotus|Lappet-faced Vulture|2|4|



##Figures

![](proposal_tree.png)
Figure 1. Hierarchical clustering of the movement characteristics of several avian species (n=8). Numbers indicate unique individuals used in the analyses. The red box indicates a clade containing all Oilbirds (n=15) and the blue box indicates a clade containing Knob-billed ducks (n=7), indicating movement has low variation and is likely stereotypical in these species. No other species clustered together, suggesting individuals of these species have greater variation in their movement characteristics. 


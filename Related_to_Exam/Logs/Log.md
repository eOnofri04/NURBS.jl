# 13-novembre-2018

Dal primo incontro è risultata una suddivisione parziale delle task.

Al momento abbiamo deciso di implementare in maniera individuale (entro giovedì 15/11) gli algoritmi:
 - `basis`
 - `dbasis`
 - `dbasisu`
 - `knot`
 - `knotu`

Ovvero i 5 algoritmi di base per la crezione dei restanti 22.

Eseguito questo lavoro individuale, il passo successivo sarà quello di confrontarci su questa implementazione
e dividere le task in base alle difficoltà e alle criticità riscontrate.

# 15-novembre-2018

Ho parlato con il professor D'Autilia che mi ha suggerito di organizzare il lavoro nella seguente maniera.
 - Lo scopo è quello di creare un package (quindi a breve, nel we, riorganizzerò la struttura delle cartelle e dei file)
 - Scrivere un file `.jl` di gestione e un file `.jl` per ogni funzione.
 - L'idea è quella di sviluppare alcuni esempi di superfici (ad esempio superfici a curvatura media costanti) e di porci quello come obbiettivo ultimo.
 - La suddivisione del lavoro non sarà ben determinata ma sarà di tipo _WorkFlow_, ovvero quando uno _**vuole**_ lavorare su una funzione:
   1. Scrive sul file di log di cosa si sta occupando (così da non entrare in conflitto con gli altri)
   1. Effettua un branching del progetto
   1. Lavora su ciò che ha dichiarato individualmente
   1. A fine lavoro ricollega il brach
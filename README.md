# DendroCrhAn

Questo repository contiene i primi risultati dell'elaborazione delle serie dendrocronologiche europee.

<mark>VISIONARE FILE esempio.pdf PER UN ESEMPIO DI CARICAMENTO E UTILIZZO DEI DATI</mark>

## Contenuto dei file
I file sono salvati in formato `.rds` per essere facilmente caricati in R.

* **`cronologie_aggiornate.rds`**: Contiene la lista delle medie dei campi delle serie dendrocronologiche detrendizzate
* **`mappa_stand_serie_storiche.rds`**: Contiene l'oggetto spaziale dei punti, un raster dove sono segnate le coordinate, e per ciascuna coppia è contenuta la serie di crescita del campo.
* **`mappa_stand_serie_storiche.rds`**: Contiene un data frame con i metadata di ciascun sito: "stand_id", "longitude", "latitude" "altitude", "zona_id", "climate_zona", "species", "header_text", le informazioni estratte (a parte latitudine, longitudine e specie) non sono sempre presenti, motivo per cui ho voluto includere tutta l'intestazione del file del sito ("header_text") per poter recuperare informazioni che non sono state catturate. Le coordinate geografiche sono invece state controllate insieme ai dati climatici e corrette nel caso risultassero inesattezze.
* **`pack.RData`**: Contiene un salvataggio di un Environment di R, simulando un pacchetto contenente tutte le funzioni scritte da me e utilizzate per svolgere le analisi (purtroppo non sono presenti descrizioni delle funzioni)

## Come caricare i file in R
Per riutilizzare questi file nel tuo script, usa il seguente codice:

```R
# Per le cronologie
risultati <- readRDS("cronologie_aggiornate.rds")

# Per la mappa (richiede il pacchetto terra)
library(terra)
mappa <- unwrap(readRDS("mappa_stand_serie_storiche.rds"))

# Per i metadati
df_meta <- readRDS("metadata_zone_climatiche.rds")

# Per il pacchetto di funzioni
load("pack.RData")

```

## Fonti dei dati
Dati dendrocronologici: https://www.ncei.noaa.gov/pub/data/paleo/treering/measurements/europe/

Zonce climatiche: https://koeppen-geiger.vu-wien.ac.at/present.htm

Dati climatici SPEI: https://spei.csic.es/database.html

Altri dati climatici: https://catalogue.ceda.ac.uk/uuid/9cf07e92afaa405da4f40b6733f362d3/




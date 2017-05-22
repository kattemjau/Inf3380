Oblig2 inf3380
Sindre Hovland
sindrech

Kjør programmet med: make

Programmet basserer seg på antagelsesen av at kvadratroten av antall traader er et heltall.
Dvs at antall prosesser er feks 4, 9, eller 16 + main prosessen som leser fra fil, sender oppgaver til barneprosessene og skriver til fil.
Totalt vil programmet kjøre med 5, 10, 17 traader osv.
programmet vil også dele inn oppgavene til hver barne prosess i 4 deler med open-mp

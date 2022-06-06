# WBO_Project_2
Project 2 from the course "Wstep do biologii obliczeniowej".
The author of the task: Bartosz Wilczyński, University of Warsaw
Task:

1. (7 pkt) Napisz program extend.py, który dla zadanej listy fragmentów białek w formacie fasta, znajdzie w bazie NCBI (protein non redundant)  przy pomocy programu blastp wszystkie sekwencje białkowe o zadanym minimalnym procencie identyczności (domyślnie 90%), i zadanym e-value nie większym niż zadana wartość (domyślnie 10e-10) z co najmniej jednym fragmentem z danych wejściowych.  Przykładowy plik z danymi: input-z2.fasta. Program powinien przyjmować parametry z linii komend, dokonywać zapytania przez internet, parsować wyjście i zwracać plik .fasta

2.   (7 pkt) Napisz program scan_pfam.py, który na podstawie pliku z białkami wykonuje zapytanie do serwera hmmscan (serwer online tu, opis api) i pobiera pliki wynikowe w formacie tsv. Korzystając z plików pobranych z serwera hmmer poda nam w wyniku plik csv, w którym będziemy mieli w wierszach kolejne identyfikatory białek z pliku FASTA, zaś w kolumnach będzie miał kolejne identyfikatory domen białkowych PFAM. Na przecięciu wiersza i kolumny stawiamy 0, jeśli dana domena nie została znaleziona w danym białku, a 1 w przeciwnym wypadku. Niestety moduły do parsowania wyjścia z HMMera w biopythonie często sprawiają problemy, dlatego polecam ręcznie wczytywać pliki tsv.

3.  (6 pkt ) Napisz program, który dla zadanych dwóch równych plików wynikowych z punktu 2, ( w naszym przykładzie możemy wyliczyć takie pliki dla pliku źródłowego input-z2.fasta i dla pliku wynikowegoz z punktu 1.) wyliczy nam wartość prawdopodobieństwa takiego rozkładu “trafień” przy założeniu, że w obu grupach białek domeny powinny występować równie często. W tym podpunkcie chodzi nam o napisanie implementacji testu Fishera dla naszych danych korzystając z wzoru:

![equation](./equation.png)

 , i porównanie wyników do implementacji testu Fishera zaimplementowanego w module Scipy.stats. (Oznaczenia: N – liczba białek w obu plikach, K -liczba białek w pierwszym pliku, n – liczba “trafień” w obu plikach, x – liczba “trafień” w pierwszym pliku).

Jako wynik przesyłamy skomentowane pliki .py

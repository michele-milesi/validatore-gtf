# Progetto Elementi di Bioinformatica - Validatore File GTF

La tabella sottostante specifica le violazioni considerate e gli ID ad esse associati.<br /><br />
Nella cartella tests sono presenti i file GTF che includono le violazioni considerate: il file associato ad una determinata violazoine è 'ID.gtf':
  per esempio il file che contiene la violazione 1 (Separatore illegale) è il file '1.gtf' nella cartella tests.<br /><br />
L'output della validazione è inserito nel file 'risultato.txt'. <br />


<table>
  <tr><th>ID</th><th>Violazione</th></tr>
  <tr><td>1</td><td>Il separatore dei vari campi della riga non è diverso da '\t' (singolo carattere di tabulazione)</td></tr>
  <tr><td>2</td><td>Il record obbligatorio 'start_codon' non è presente nel file</td></tr>
  <tr><td>3</td><td>Il record obbligatorio 'stop_codon' non è presente nel file</td></tr>
  <tr><td>4</td><td>Il campo source non è unico nel file</td></tr>
  <tr><td>5</td><td>'start_codon' e 'stop_codon' sono lunghi più di 3 bp</td></tr>
  <tr><td>6</td><td>'start_codon' non è incluso nelle coordinate di nessun record 'CDS'</td></tr>
  <tr><td>7</td><td>'stop_codon' è incluso nelle coordinate di almeno un record 'CDS' o '3UTR'</td></tr>
  <tr><td>8</td><td>Il numero di campi è errato: devono esserci esattamente 9 campi separati da '\t'</td></tr>
  <tr><td>9</td><td>Il separatore degli attributi di un regord è sbagliato: ogni attributo deve terminare con punto e virgola e deve essere separato dal successivo esattamente da uno spazio</td></tr>
  <tr><td>10</td><td>Attributo ha formato sbagliato: ogni attributo deve essere una coppia nome-valore, gli elementi devone essere separati esattamente da uno spazio</td></tr>
  <tr><td>11</td><td>Valore testuale di un attributo non è compreso da doppi apici</td></tr>
  <tr><td>12</td><td>Ordine attributi sbagliato: transfert_id e gene_id devono essere i primi due attributi per ogni record</td></tr>
  <tr><td>13</td><td>L'attributo obbligatorio gene_id non è presente in un record</td></tr>
  <tr><td>14</td><td>L'attributo obbligatorio transfert_id non è presente in un record</td></tr>
  <tr><td>15</td><td>Warning: Il nome degli attributi non deve essere compreso in doppi apici, ciò è dovuto al fatto che se un attributo obbligatrio fosse compreso in doppi apici, esso sarebbe visto dal validatore come un attributo diverso, ciò porebbe provocare altre violazioni</td></tr>
  <tr><td>16</td><td>'start_codon' e 'stop_codon' contigui non hanno frame uguale a 0: per definizione se sono presenti 'start_codon' o 'stop_codon' contigui allora il valore del loro frame deve essere 0</td></tr>
  <tr><td>17</td><td>'start_codon' e 'stop_codon' non contigui hanno valore del frame errato: se non sono contigui il valore del frame può essere 0 o 1 o 2</td></tr>
  <tr><td>18</td><td>Il valre del campo frame per record che non sono 'start_codon' e 'stop_codon' è errato: il valore del frame può essere '.' o 0 o 1 o 2</td></tr>
  <tr><td>19</td><td>Il valore del campo strand è errato: deve essere o '+' o '-'</td></tr>
  <tr><td>20</td><td>Il valore del campo score è illegale: può essere '.' o un intero o un floating point</td></tr>
  <tr><td>21</td><td>Il valore del campo start è illegale: deve essere un intero maggiore o uguale a 1</td></tr>
  <tr><td>22</td><td>Il valore del campo end è illegale: deve essere un intero maggiore o uguale a 1</td></tr>
  <tr><td>23</td><td>Il valore del vampo start è maggiore del valore del campo end</td></tr>
  <tr><td>24</td><td>'start_codon' incluso nelle coordinate di un record '5UTR'</td></tr>
  <tr><td>25</td><td>L'attributo transcript_id per il record 'intron_CNS' è vuoto</td></tr>
  <tr><td>26</td><td>L'attributo transcript_id per i record 'inter' e 'inter_CNS' non è vuoto</td></tr>
  <tr><td>27</td><td>Il record obbligatorio 'CDS' non è presente nel file</td></tr>
</table>

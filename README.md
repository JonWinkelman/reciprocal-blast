## Reciprocal blast

<img src="images/reciprocal_blast.png" alt="reciprocal blast">

dependencies:<br>
python3
`pandas` `Bio` `plotly`<br> 
BLAST+  <a href="https://www.ncbi.nlm.nih.gov/books/NBK569861/#intro_Installation.MacOSX">Install</a>
<br>
<a href="https://www.ncbi.nlm.nih.gov/books/NBK279691/">BLAST+ user manual</a>
<br>
1) build individual local blast database from each of the proteomes<br>  
2) blast each proteome database with an input fasta query protein<br>   
3) turn best hit from each result in 2) into an individual single best hit fasta <br>  
4) run each best hit from 3) against the genome containing the original query sequence<br>  
5) build csv for results in 4<br> 
<br>

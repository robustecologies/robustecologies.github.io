# Detailed diagram of the Samkhya derivation

Convenience wrapper for `create_samkhya_flowchart(detailed = TRUE)`,
which labels every node with its index range and names each principle it
contains.

## Usage

``` r
create_samkhya_flowchart_detailed()
```

## Value

A character string containing a Mermaid `graph TD` definition.

## Details

The numbering is that of
[tattvas](https://robustecologies.github.io/SamkhyaR/reference/tattvas.md):
the twenty-four material principles run from mulaprakrti (1) to prthivi
(24) and Purusha is the twenty-fifth. Kaivalya carries no number,
because it is a condition of Prakrti and of buddhi and not a principle
of the enumeration.

## References

Isvarakrsna. *Samkhyakarika*, verses 3, 22, 25-26, 38.

## See also

[`create_samkhya_flowchart()`](https://robustecologies.github.io/SamkhyaR/reference/create_samkhya_flowchart.md)
of which this is the detailed form,
[`print_samkhya_flowchart()`](https://robustecologies.github.io/SamkhyaR/reference/print_samkhya_flowchart.md),
[tattvas](https://robustecologies.github.io/SamkhyaR/reference/tattvas.md),
[`plot.samkhya_state()`](https://robustecologies.github.io/SamkhyaR/reference/plot.samkhya_state.md).

## Examples

``` r
cat(create_samkhya_flowchart_detailed())
#> graph TD
#>     P(("<b>25. Purusha</b><br/>nirguna, akartr, saksin"))
#>     A["<b>1. Mulaprakriti</b><br/>avyakta, the uncreated root"]
#>     B["<b>2. Mahat / Buddhi</b><br/>adhyavasaya, ascertainment"]
#>     C{"<b>3. Ahamkara</b><br/>abhimana, appropriation"}
#>     T["<b>taijasa</b><br/>rajasic aspect<br/>emits no tattva"]
#> 
#>     subgraph vaikrta aspect - sattvic - the eleven organs
#>     D["<b>4. Manas</b><br/>samkalpaka, and eleventh organ"]
#>     E["<b>5-9. Jnanendriyas</b><br/>Shrotra, Tvac, Chakshus<br/>Rasana, Ghrana"]
#>     F["<b>10-14. Karmendriyas</b><br/>Vac, Pani, Pada<br/>Payu, Upastha"]
#>     end
#> 
#>     subgraph bhutadi aspect - tamasic - the elements
#>     G["<b>15-19. Tanmatras</b><br/>Shabda, Sparsha, Rupa<br/>Rasa, Gandha<br/>avisesa, the non-specific"]
#>     H["<b>20-24. Mahabhutas</b><br/>Akasha, Vayu, Tejas<br/>Ap, Prithivi<br/>visesa, the specific"]
#>     end
#> 
#>     K["<b>Kaivalya</b><br/>Prakrti withdraws, having been seen<br/>SK 59, 66-68. A condition, not a tattva."]
#> 
#>     P -. "sannidhimatra: presence, not action. SK 21" .-> A
#>     A --> |SK 22| B
#>     B --> |SK 22| C
#>     C --> |vaikrta. SK 25| D
#>     C --> |vaikrta. SK 25| E
#>     C --> |vaikrta. SK 25| F
#>     C --> |bhutadi. SK 25| G
#>     G --> |SK 22, 38| H
#>     T -. drives .-> D
#>     T -. drives .-> G
#>     H -. "viveka arises in buddhi. SK 64" .-> K
#> 
#>     style P fill:#FFFFFF,stroke:#170C3A,stroke-width:3px,color:#000
#>     style A fill:#4ECDC4,stroke:#0D7377,stroke-width:2px,color:#000
#>     style B fill:#E4F1E8,stroke:#38A169,stroke-width:2px,color:#000
#>     style C fill:#F0D9DC,stroke:#A6081A,stroke-width:3px,color:#000
#>     style T fill:#A6081A,stroke:#6A0410,stroke-width:2px,color:#fff
#>     style D fill:#F2F2F2,stroke:#5A3FA0,stroke-width:2px,color:#000
#>     style E fill:#F2F2F2,stroke:#5A3FA0,stroke-width:2px,color:#000
#>     style F fill:#F2F2F2,stroke:#5A3FA0,stroke-width:2px,color:#000
#>     style G fill:#4A4066,stroke:#170C3A,stroke-width:2px,color:#fff
#>     style H fill:#170C3A,stroke:#000000,stroke-width:2px,color:#fff
#>     style K fill:#0054AD,stroke:#00306B,stroke-width:3px,color:#fff
```

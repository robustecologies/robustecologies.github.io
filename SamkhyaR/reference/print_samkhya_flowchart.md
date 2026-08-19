# Render the Samkhya diagram

Returns the Mermaid definition of the derivation diagram and, when
`DiagrammeR` is available and rendering is requested, draws it.

## Usage

``` r
print_samkhya_flowchart(
  detailed = FALSE,
  render = interactive(),
  verbose = FALSE
)
```

## Arguments

- detailed:

  Logical; use the detailed labelling. Default `FALSE`.

- render:

  Logical; draw the diagram with
  [`DiagrammeR::mermaid()`](https://rich-iannone.github.io/DiagrammeR/reference/mermaid.html).
  Default `TRUE` in an interactive session and `FALSE` otherwise, so
  that the function is silent in scripts and in checks.

- verbose:

  Logical; print the Mermaid source and a note on how to render it
  elsewhere. Default `FALSE`.

## Value

The Mermaid definition, invisibly.

## Details

`DiagrammeR` is required only for `render = TRUE`; the definition itself
is a plain character string that any Mermaid renderer will accept,
including a `mermaid` chunk in an R Markdown document and the online
editor at <https://mermaid.live>.

## References

Isvarakrsna. *Samkhyakarika*, verses 22, 25.

## See also

[`create_samkhya_flowchart()`](https://robustecologies.github.io/SamkhyaR/reference/create_samkhya_flowchart.md)
and
[`create_samkhya_flowchart_detailed()`](https://robustecologies.github.io/SamkhyaR/reference/create_samkhya_flowchart_detailed.md)
for the definitions themselves,
[`plot.samkhya_state()`](https://robustecologies.github.io/SamkhyaR/reference/plot.samkhya_state.md)
for the ggplot2 rendering.

## Examples

``` r
definition <- print_samkhya_flowchart(render = FALSE)
substr(definition, 1, 8)
#> [1] "graph TD"

print_samkhya_flowchart(detailed = TRUE, render = FALSE, verbose = TRUE)
#> Samkhya derivation diagram (detailed labelling):
#> 
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
#> 
#> To render elsewhere, paste into a mermaid chunk or into https://mermaid.live
```

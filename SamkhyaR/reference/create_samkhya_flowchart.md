# Diagram of the Samkhya derivation

Builds a Mermaid graph of the derivation of the twenty-five principles,
following Samkhya-Karika 22 and 25. The graph is generated from
[tattvas](https://robustecologies.github.io/SamkhyaR/reference/tattvas.md),
so its numbering and its structure cannot drift from the enumeration the
rest of the package uses.

## Usage

``` r
create_samkhya_flowchart(detailed = FALSE)
```

## Arguments

- detailed:

  Logical; name every member of each group and give its index range,
  rather than labelling the group alone. Default `FALSE`.

## Value

A character string containing a Mermaid `graph TD` definition, suitable
for
[`DiagrammeR::mermaid()`](https://rich-iannone.github.io/DiagrammeR/reference/mermaid.html)
or for a `mermaid` chunk in an R Markdown document.

## Details

The two aspects of ahamkara are drawn as two enclosed branches, because
their parallelism is the point of SK 25 and the point most often lost in
diagrams of this system. The *vaikrta* branch carries manas and the ten
organs; the *bhutadi* branch carries the five tanmatras, from which
alone the five gross elements descend. The *taijasa* aspect is drawn as
a node attached to both by broken edges, since it supplies the activity
by which either proceeds and emits no principle of its own.

Node colour encodes the quality that predominates: pale for the sattvic
branch, dark for the tamasic, red for the rajasic aspect that drives
both. Purusha is drawn uncoloured, outside the sequence, and joined to
it by a broken edge, because it is *nirguna*, constituted by none of the
three, and because SK 3 places it outside the causal series as neither
cause nor effect. Its relation to Prakrti is *sannidhimatra*, mere
proximity: it does not act (SK 21).

Three things this diagram deliberately does not do. It does not route
the organs into the tanmatras, which would invert SK 25. It does not
assign the organs of action to a tamasic or rajasic branch, which would
contradict the same verse. And it does not present kaivalya as a
twenty-fifth principle: the twenty-fifth principle is Purusha, and
kaivalya is a condition of Prakrti and of buddhi, not a tattva.

## References

Isvarakrsna. *Samkhyakarika*, verses 3, 21-22, 25-26, 38, 59, 64, 66-68.

## See also

[`print_samkhya_flowchart()`](https://robustecologies.github.io/SamkhyaR/reference/print_samkhya_flowchart.md)
for a wrapper that renders it,
[`plot.samkhya_state()`](https://robustecologies.github.io/SamkhyaR/reference/plot.samkhya_state.md)
for a state-aware rendering of the same graph,
[tattvas](https://robustecologies.github.io/SamkhyaR/reference/tattvas.md)
for the enumeration it is generated from,
[`samkhya_evolve()`](https://robustecologies.github.io/SamkhyaR/reference/samkhya_evolve.md).

## Examples

``` r
cat(create_samkhya_flowchart())
#> graph TD
#>     P(("<b>Purusha</b><br/>nirguna, akartr, saksin"))
#>     A["<b>Mulaprakriti</b><br/>avyakta, the uncreated root"]
#>     B["<b>Mahat / Buddhi</b><br/>adhyavasaya, ascertainment"]
#>     C{"<b>Ahamkara</b><br/>abhimana, appropriation"}
#>     T["<b>taijasa</b><br/>rajasic aspect<br/>emits no tattva"]
#> 
#>     subgraph vaikrta aspect - sattvic - the eleven organs
#>     D["<b>Manas</b><br/>samkalpaka, and eleventh organ"]
#>     E["<b>Jnanendriyas</b><br/>five organs of knowledge<br/>alocanamatra, bare apprehension"]
#>     F["<b>Karmendriyas</b><br/>five organs of action"]
#>     end
#> 
#>     subgraph bhutadi aspect - tamasic - the elements
#>     G["<b>Tanmatras</b><br/>five bare essences<br/>avisesa, the non-specific"]
#>     H["<b>Mahabhutas</b><br/>five gross elements<br/>visesa, the specific"]
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

# The detailed form names every principle.
cat(create_samkhya_flowchart(detailed = TRUE))
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

# \donttest{
if (requireNamespace("DiagrammeR", quietly = TRUE)) {
    DiagrammeR::mermaid(create_samkhya_flowchart())
}

{"x":{"diagram":"graph TD\n    P((\"<b>Purusha<\/b><br/>nirguna, akartr, saksin\"))\n    A[\"<b>Mulaprakriti<\/b><br/>avyakta, the uncreated root\"]\n    B[\"<b>Mahat / Buddhi<\/b><br/>adhyavasaya, ascertainment\"]\n    C{\"<b>Ahamkara<\/b><br/>abhimana, appropriation\"}\n    T[\"<b>taijasa<\/b><br/>rajasic aspect<br/>emits no tattva\"]\n\n    subgraph vaikrta aspect - sattvic - the eleven organs\n    D[\"<b>Manas<\/b><br/>samkalpaka, and eleventh organ\"]\n    E[\"<b>Jnanendriyas<\/b><br/>five organs of knowledge<br/>alocanamatra, bare apprehension\"]\n    F[\"<b>Karmendriyas<\/b><br/>five organs of action\"]\n    end\n\n    subgraph bhutadi aspect - tamasic - the elements\n    G[\"<b>Tanmatras<\/b><br/>five bare essences<br/>avisesa, the non-specific\"]\n    H[\"<b>Mahabhutas<\/b><br/>five gross elements<br/>visesa, the specific\"]\n    end\n\n    K[\"<b>Kaivalya<\/b><br/>Prakrti withdraws, having been seen<br/>SK 59, 66-68. A condition, not a tattva.\"]\n\n    P -. \"sannidhimatra: presence, not action. SK 21\" .-> A\n    A --> |SK 22| B\n    B --> |SK 22| C\n    C --> |vaikrta. SK 25| D\n    C --> |vaikrta. SK 25| E\n    C --> |vaikrta. SK 25| F\n    C --> |bhutadi. SK 25| G\n    G --> |SK 22, 38| H\n    T -. drives .-> D\n    T -. drives .-> G\n    H -. \"viveka arises in buddhi. SK 64\" .-> K\n\n    style P fill:#FFFFFF,stroke:#170C3A,stroke-width:3px,color:#000\n    style A fill:#4ECDC4,stroke:#0D7377,stroke-width:2px,color:#000\n    style B fill:#E4F1E8,stroke:#38A169,stroke-width:2px,color:#000\n    style C fill:#F0D9DC,stroke:#A6081A,stroke-width:3px,color:#000\n    style T fill:#A6081A,stroke:#6A0410,stroke-width:2px,color:#fff\n    style D fill:#F2F2F2,stroke:#5A3FA0,stroke-width:2px,color:#000\n    style E fill:#F2F2F2,stroke:#5A3FA0,stroke-width:2px,color:#000\n    style F fill:#F2F2F2,stroke:#5A3FA0,stroke-width:2px,color:#000\n    style G fill:#4A4066,stroke:#170C3A,stroke-width:2px,color:#fff\n    style H fill:#170C3A,stroke:#000000,stroke-width:2px,color:#fff\n    style K fill:#0054AD,stroke:#00306B,stroke-width:3px,color:#fff\n"},"evals":[],"jsHooks":[]}# }
```

# Changelog

## SamkhyaR 0.2.0

### Breaking changes

- The evolution functions now return an S3 object of class
  `samkhya_state` rather than a bare list. The field `tattvas_evolved`
  is replaced by `tattvas`, which holds IAST names indexed into the new
  `tattvas` data object, and `gunas` is now a named numeric vector
  rather than a list. Code that reached into the old list will need
  updating.

- [`evolution_buddhi()`](https://robustecologies.github.io/SamkhyaR/reference/evolution_buddhi.md)
  validates `guna_dominant` with
  [`match.arg()`](https://rdrr.io/r/base/match.arg.html). An
  unrecognised value previously passed silently and had no effect; it is
  now an error.

- [`print_samkhya_flowchart()`](https://robustecologies.github.io/SamkhyaR/reference/print_samkhya_flowchart.md)
  gains `render` and `verbose` arguments, both of which default to not
  writing to the console outside an interactive session. It previously
  printed unconditionally.

### Bug fixes

- [`get_kaivalya()`](https://robustecologies.github.io/SamkhyaR/reference/get_kaivalya.md)
  no longer returns its argument unchanged. The function gated on the
  evolved principles numbering 15 while the documented pipeline produced
  8, so the condition could never be satisfied: liberation silently did
  nothing, no message was emitted, and `final_state` was `NULL`. The
  function now gates on the twenty-three evolutes of Prakrti actually
  being present, and errors with a diagnostic when discrimination is
  incomplete.

- The three qualities now remain on the two-simplex.
  [`evolution_buddhi()`](https://robustecologies.github.io/SamkhyaR/reference/evolution_buddhi.md)
  added 0.1 to the dominant quality without removing it from the others,
  so the triple summed to 1.1 after the first step and drifted further
  with each disturbance. The disturbance is now a redistribution, and
  every function validates the invariant.

- The derivation diagram no longer contradicts Samkhya-Karika 25. The
  organs of action were labelled a tamasic branch in one diagram and a
  rajasic branch in the other, where the verse assigns all eleven organs
  to the sattvic (vaikrta) aspect of ahamkara; and both diagrams routed
  manas and the organs into the tanmatras, where the verse derives the
  tanmatras from the tamasic (bhutadi) aspect of ahamkara in parallel
  with the organs. Both diagrams are now generated from the `tattvas`
  data object.

- [`evolution_tanmatras()`](https://robustecologies.github.io/SamkhyaR/reference/evolution_tanmatras.md)
  no longer presupposes that the organs have manifested. The two
  branches from ahamkara are parallel, so a state carrying ahamkara
  alone is a valid argument.

- The count of Prakrti’s evolutes is corrected from 24 to 23 throughout.
  Mulaprakrti is not an evolute of herself (SK 3).

- The subtle body is no longer said to be constituted by the gross
  elements. SK 40 makes it the eighteen consisting of mahat, ahamkara,
  the eleven organs and the five tanmatras.

- Verse attributions are corrected where the karika does not say what
  was claimed: SK 27 gives `samkalpaka` and not `savikalpa`; SK 23
  reverses `viraga` to `raga` and not to `avairagya`; SK 38 predicates
  tranquil, turbulent and deluding of the gross elements and not of the
  tanmatras; and the cumulative-quality scheme of the gross elements is
  Vijnanabhiksu’s position rather than the text of SK 22.

- The tattva numbering no longer differs between the reference pages and
  the README, and the organs of knowledge are listed in one order
  throughout rather than three.

- [`get_kaivalya()`](https://robustecologies.github.io/SamkhyaR/reference/get_kaivalya.md)
  no longer alters the description of Purusha. SK 62 denies that anyone
  is bound, released or a transmigrant, so what the function records is
  the arising of discriminative knowledge in buddhi and the cessation of
  Prakrti’s activity.

### New features

- New data object `tattvas`, the twenty-five principles with the
  structural partition of SK 3, the immediate material cause of each,
  the aspect of ahamkara from which it issues under SK 25, and its
  membership in the internal organ, the elevenfold set and the subtle
  body.

- New S3 class `samkhya_state` with
  [`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rdrr.io/r/base/summary.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods. The
  plot method offers `type = "gunas"`, `"tattvas"` and `"hierarchy"`.

- New function
  [`samkhya_evolve()`](https://robustecologies.github.io/SamkhyaR/reference/samkhya_evolve.md)
  runs the complete branching sequence in one call.

- New doctrinal reference functions, each returning a classed object
  with a [`print()`](https://rdrr.io/r/base/print.html) method:
  [`samkhya_pramanas()`](https://robustecologies.github.io/SamkhyaR/reference/samkhya_pramanas.md)
  for the three means of valid knowledge in the karika’s own vocabulary
  (SK 4-6),
  [`buddhi_bhavas()`](https://robustecologies.github.io/SamkhyaR/reference/buddhi_bhavas.md)
  for the eight dispositions and the fiftyfold intellectual creation (SK
  23, 43-51),
  [`linga_sharira()`](https://robustecologies.github.io/SamkhyaR/reference/linga_sharira.md)
  for the eighteen constituents of the subtle body (SK 39-42), and
  [`satkaryavada()`](https://robustecologies.github.io/SamkhyaR/reference/satkaryavada.md)
  for the five arguments of SK 9.

- [`evolution_indriyas()`](https://robustecologies.github.io/SamkhyaR/reference/evolution_indriyas.md)
  gains a `which` argument, so the organs of knowledge and of action can
  be manifested separately.

- [`get_kaivalya()`](https://robustecologies.github.io/SamkhyaR/reference/get_kaivalya.md)
  gains a `jivanmukti` argument distinguishing the persistence of the
  body by momentum (SK 67) from the final separation (SK 68).

- Every evolution function gains a `verbose` argument.

- [`create_samkhya_flowchart()`](https://robustecologies.github.io/SamkhyaR/reference/create_samkhya_flowchart.md)
  gains a `detailed` argument;
  [`create_samkhya_flowchart_detailed()`](https://robustecologies.github.io/SamkhyaR/reference/create_samkhya_flowchart_detailed.md)
  is retained as a wrapper.

### Documentation

- Examples now run rather than sitting inside `\dontrun{}`, so
  `R CMD check` exercises them.

## SamkhyaR 0.1.0

- Initial release.

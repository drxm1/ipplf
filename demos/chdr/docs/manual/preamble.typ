// One import provides the manual helpers and all physica notation.
#import "@preview/physica:0.9.8": *

// Shared typography for the full manual and standalone chapters.
#let manual-typography(body) = {
  set text(font: "Libertinus Serif", size: 11pt, lang: "en", fill: rgb("202428"))
  set par(justify: true, leading: 0.6em)
  set heading(numbering: "1.1", supplement: [Section])
  set math.equation(numbering: "(1)")
  show math.equation: set text(font: "New Computer Modern Math")
  show heading: set text(weight: "semibold")
  show heading.where(level: 1): set text(size: 16pt)
  show heading.where(level: 2): set text(size: 12pt)
  show heading.where(level: 3): set text(size: 11pt)
  show raw: set text(font: "DejaVu Sans Mono", size: 0.85em)
  body
}

// Keep document metadata and page geometry out of chapter bodies.
#let manual-layout(body, title: "ChDR Manual") = {
  set document(title: title)
  set page(
    paper: "a4",
    margin: (x: 24mm, y: 24mm),
    footer: context align(center, text(
      size: 9pt,
      fill: luma(45%),
      counter(page).display("1"),
    )),
  )
  manual-typography(body)
}

#let manual-title(title) = block(above: 0pt, below: 18pt, stack(
  dir: ttb,
  spacing: 8pt,
  text(size: 24pt, weight: "semibold", title),
  line(length: 100%, stroke: 0.45pt + luma(80%)),
))

// Always register the bibliography so @key references can resolve. Only its
// heading is conditional; an unused bibliography contributes no visible text.
#let manual-references() = context {
  bibliography(
    "references.bib",
    style: "ieee",
    title: if query(cite).len() > 0 [References] else { none },
  )
}

// The root controls chapter order. Enable the outline once it is useful.
#let manual(body, title: "ChDR Manual", contents: false) = manual-layout(
  title: title,
)[
  #manual-title(title)
  #if contents { outline(title: [Contents], depth: 2) }
  #body
  #manual-references()
]

// A chapter supplies its own headings. Its standalone rendering is discarded
// when the root imports only the exported chapter content.
#let standalone(body, title: "ChDR Manual") = manual-layout(
  title: title,
)[
  #body
  #manual-references()
]

// One import provides the manual helpers and all physica notation.
#import "@preview/physica:0.9.8": *
#import "@preview/equate:0.3.3": equate
#import "@preview/simple-todo:0.1.0": list-todos, todo

#let proposal(body) = text(fill: rgb("1f4e9c"), body)
#let revision(body) = text(fill: rgb("707070"), body)

// Shared typography for the manual.
#let manual-typography(body) = {
  set text(
    font: "Libertinus Serif",
    size: 11pt,
    lang: "en",
    fill: rgb("202428"),
  )
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
    width: 210mm,
    height: auto,
    margin: (x: 24mm, y: 24mm),
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
  // Temporary review markup for the added bibliography entries.
  let proposed = cite.where(key: <xxxx>).or(cite.where(key: <yyyy>))
  bibliography(
    "references.bib",
    style: "ieee",
    title: if query(cite).len() > 0 [References] else { none },
  )
  if query(proposed).len() > 0 {
    block(
      width: 100%,
      inset: 10pt,
      stroke: 0.7pt + rgb("315f86"),
      fill: rgb("f4f8fc"),
    )[
      #set text(fill: rgb("234e70"))
      #bibliography(
        "references.bib",
        style: "ieee",
        title: none,
        target: proposed,
      )
    ]
  }
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

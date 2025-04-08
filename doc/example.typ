#import "lib.typ": *
#show: template.with(
  title: [TightBinding Model Document],
  short_title: "manual",
  description: [
    A manual for the Slater-Koster TightBinding Model
  ],
  authors: (
    (
      name: "Quinn Hsu",
      link: "https://github.com/HsuQuinn",
    ),
  ),
  paper_size: "a4",
  // landscape: true,
  cols: 1,
  
  text_font: "Proxima Nova",
  code_font: "Cascadia Mono",
  accent: "#15AD66", // green
  h1-prefix: "Section",
  colortab: true,
)

#include "content/doc.typ"

// Import type definitions to make development easier
// These import lines are commented out at bundling into a single HTML file
import MicroModal from "./lib/types/micromodal"
import { RowComponent, Tabulator } from "./lib/types/tabulator"

const FEATURES_JSON = {}
const LINKS_JSON = {}

const SVG_ANIMATE_HTML =
  "<animate attributeName='fill-opacity', begin='0s', dur='1s' values='0.3;1;0.3' repeatCount='indefinite' />"

/**
 * Download PNG image
 * @param {HTMLElement} svgNode
 * @param {string} fileName
 */
function downloadPng(svgNode, fileName = "image.png") {
  const svgData = new XMLSerializer().serializeToString(svgNode)
  const canvas = document.createElement("canvas")
  canvas.width = svgNode.width.baseVal.value
  canvas.height = svgNode.height.baseVal.value

  const ctx = canvas.getContext("2d")
  const image = new Image()
  image.onload = function () {
    ctx.drawImage(image, 0, 0)
    const a = document.createElement("a")
    a.href = canvas.toDataURL("image/png")
    a.setAttribute("download", fileName)
    a.dispatchEvent(new MouseEvent("click"))
  }
  image.src =
    "data:image/svg+xml;charset=utf-8;base64," +
    btoa(decodeURIComponent(encodeURIComponent(svgData)))
}
/**
 * Download SVG image
 * @param {HTMLElement} svgNode
 * @param {string} fileName
 */
function downloadSvg(svgNode, fileName = "image.svg") {
  const svgData = new XMLSerializer().serializeToString(svgNode)
  const svgBlob = new Blob([svgData], { type: "image/svg+xml;charset=utf-8" })

  const a = document.createElement("a")
  a.href = URL.createObjectURL(svgBlob)
  a.setAttribute("download", fileName)
  a.dispatchEvent(new MouseEvent("click"))
}

/**
 * Get gid & feature path list object from SVG path list
 * (Exon-Intron feature may contain multiple path)
 * @param {HTMLCollectionOf<SVGPathElement>} svgPathList
 * @returns {Object.<string, SVGPathElement[]} gid & featurePathList
 */
function getGid2FeaturePathList(svgPathList) {
  let gid2FeaturePathList = {}
  for (let path of svgPathList) {
    const gid = path.parentNode.id
    if (gid.startsWith("Feature")) {
      if (gid in gid2FeaturePathList) {
        gid2FeaturePathList[gid].push(path)
      } else {
        gid2FeaturePathList[gid] = [path]
      }
    }
  }
  return gid2FeaturePathList
}

/**
 * Convert feature json to HTML table for tooltip display
 * @param {Object} featureJson
 * @returns {string} Feature HTML table
 */
function convertFeatureJsonToHtmlTable(featureJson) {
  const segment = featureJson.segment
  const location = featureJson.location
  const length = featureJson.length.toLocaleString("en-US")
  const type = featureJson.type
  const gene = featureJson.gene
  const protein_id = featureJson.protein_id
  const product = featureJson.product
  const pseudo = featureJson.pseudo
  const extra = featureJson.extra

  if (type === "na" && gene === "na" && protein_id == "na" && product == "na") {
    let tableLines = [
      // Simple Feature Table HTML
      `<table>`,
      `<tr><td><b>segment </b></td><td>${segment}</td></tr>`,
      `<tr><td><b>location </b></td><td>${location}</td></tr>`,
      `<tr><td><b>length </b></td><td>${length}</td></tr>`,
    ]
    for (let key in extra) {
      tableLines.push(`<tr><td><b>${key} </b></td><td>${extra[key]}</td></tr>`)
    }
    tableLines.push(`</table>`)
    return tableLines.join("\n")
  } else {
    let tableLines = [
      // Full Feature Table HTML
      `<table>`,
      `<tr><td><b>segment </b></td><td>${segment}</td></tr>`,
      `<tr><td><b>location </b></td><td>${location}</td></tr>`,
      `<tr><td><b>length </b></td><td>${length}</td></tr>`,
      `<tr><td><b>type </b></td><td>${type}</td></tr>`,
      `<tr><td><b>gene </b></td><td>${gene}</td></tr>`,
      `<tr><td><b>protein_id </b></td><td>${protein_id}</td></tr>`,
      `<tr><td><b>product </b></td><td>${product}</td></tr>`,
      `<tr><td><b>pseudo </b></td><td>${pseudo}</td></tr>`,
    ]
    for (let key in extra) {
      tableLines.push(`<tr><td><b>${key} </b></td><td>${extra[key]}</td></tr>`)
    }
    tableLines.push(`</table>`)
    return tableLines.join("\n")
  }
}

/**
 * Convert link json to HTML table for tooltip display
 * @param {Object} linkJson
 * @returns {string} Link HTML table
 */
function convertLinkJsonToHtmlTable(linkJson) {
  const start1 = linkJson.start1.toLocaleString("en-US")
  const start2 = linkJson.start2.toLocaleString("en-US")
  const end1 = linkJson.end1.toLocaleString("en-US")
  const end2 = linkJson.end2.toLocaleString("en-US")
  const length1 = linkJson.length1.toLocaleString("en-US")
  const length2 = linkJson.length2.toLocaleString("en-US")
  const identity = linkJson.identity === "na" ? "na" : `${linkJson.identity.toFixed(2)}%`
  return [
    // Link Table HTML
    `<table>`,
    `<tr><td><b>track </b></td><td>${linkJson.track1}</td><td>${linkJson.track2}</td></tr>`,
    `<tr><td><b>segment </b></td><td>${linkJson.segment1}</td><td>${linkJson.segment2}</td></tr>`,
    `<tr><td><b>start </b></td><td>${start1}</td><td>${start2}</td></tr>`,
    `<tr><td><b>end </b></td><td>${end1}</td><td>${end2}</td></tr>`,
    `<tr><td><b>length </b></td><td>${length1}</td><td>${length2}</td></tr>`,
    `<tr><td><b>identity </b></td><td colspan="2">${identity}</td></tr>`,
    `</table>`,
  ].join("\n")
}

/**
 * Convert tabulator feature row to fasta
 * @param {RowComponent} row
 * @returns {string} Fasta
 */
function convertFeatureRowToFasta(row) {
  let proteinId = row.getData().protein_id
  if (proteinId === "na") {
    proteinId = row.getData().track
  }
  const location = row.getData().location
  const product = row.getData().product
  const translation = row.getData().translation
  const fasta = `>${proteinId} location:${location} product:${product}\n${translation}`
  return fasta
}

/**
 * Get filter menu list for tabulator contextmenu
 * @returns {Object.<string, any>[]} Filter Menu List
 */
function getFilterMenuList() {
  let filterMenuList = []
  const targetFieldList = ["track", "segment", "type", "gene", "product"]
  for (let targetField of targetFieldList) {
    filterMenuList.push({
      label: (row) => {
        let fieldValue = String(row.getData()[targetField])
        const maxLen = 30
        if (fieldValue.length > maxLen) {
          fieldValue = fieldValue.substring(0, maxLen) + "..."
        }
        return `Filter by ${targetField}='${fieldValue}'`
      },
      action: (e, row) => {
        row.getTable().setHeaderFilterValue(targetField, row.getData()[targetField])
      },
    })
  }
  return filterMenuList
}

document.addEventListener("DOMContentLoaded", function () {
  // Initialize Modal
  MicroModal.init()

  // Setup colorpicker
  const initColor = "#ff0000" // = red
  document.getElementById("colorpicker").style.color = initColor
  document.getElementById("apply_color").style.backgroundColor = `${initColor}88`
  $("#colorpicker").spectrum({
    color: initColor,
    showInput: true,
    showPalette: true,
    showAlpha: true,
    showButtons: true,
    preferredFormat: "hex",
    hideAfterPaletteSelect: true,
    palette: [
      ["#f00", "#f90", "#ff0", "#0f0", "#0ff", "#00f", "#90f", "#f0f"],
      ["#f4cccc", "#fce5cd", "#fff2cc", "#d9ead3", "#d0e0e3", "#cfe2f3", "#d9d2e9", "#ead1dc"],
      ["#ea9999", "#f9cb9c", "#ffe599", "#b6d7a8", "#a2c4c9", "#9fc5e8", "#b4a7d6", "#d5a6bd"],
      ["#e06666", "#f6b26b", "#ffd966", "#93c47d", "#76a5af", "#6fa8dc", "#8e7cc3", "#c27ba0"],
      ["#c00", "#e69138", "#f1c232", "#6aa84f", "#45818e", "#3d85c6", "#674ea7", "#a64d79"],
      ["#900", "#b45f06", "#bf9000", "#38761d", "#134f5c", "#0b5394", "#351c75", "#741b47"],
      ["#600", "#783f04", "#7f6000", "#274e13", "#0c343d", "#073763", "#20124d", "#4c1130"],
      ["#000", "#444", "#666", "#999", "#ccc", "#eee", "#f3f3f3", "#fff"],
    ],
    change: () => {
      let pickColor = $("#colorpicker").spectrum("get").toHexString()
      document.getElementById("apply_color").style.backgroundColor = `${pickColor}88`
      if (pickColor === "#ffffff") {
        pickColor = "#ffffff00"
      }
      document.getElementById("colorpicker").style.color = pickColor
    },
  })

  const svg = document.querySelector("#svg_canvas>svg")
  const svgPathList = svg.getElementsByTagName("path")
  for (let path of svgPathList) {
    // Setup feature/link path tooltip
    const gid = path.parentNode.id
    let htmlTable = ""
    if (gid.startsWith("Feature")) {
      htmlTable = convertFeatureJsonToHtmlTable(FEATURES_JSON[gid])
    } else if (gid.startsWith("Link")) {
      htmlTable = convertLinkJsonToHtmlTable(LINKS_JSON[gid])
    } else {
      continue
    }
    tippy(`#${gid}`, {
      theme: "pygenomeviz",
      content: htmlTable,
      allowHTML: true,
      followCursor: true,
      placement: "bottom-start",
      arrow: false,
      duration: [0, 0],
      offset: [0, 15],
    })

    // Apply picked color to feature/link path
    path.addEventListener("dblclick", () => {
      path.style.fill = document.getElementById("colorpicker").style.color
    })

    // Highlight selected feature/link path
    path.addEventListener("mouseover", () => {
      path.insertAdjacentHTML("beforeend", SVG_ANIMATE_HTML)
    })
    path.addEventListener("mouseout", () => {
      path.querySelector("animate").remove()
    })
  }

  const allTexts = svg.getElementsByTagName("text")
  for (let text of allTexts) {
    // Setup text edit dialog
    text.addEventListener("contextmenu", (e) => {
      e.preventDefault()
      const targetText = e.currentTarget
      const originalText = targetText.textContent
      const originalFontSize = parseInt(targetText.style.fontSize)

      // Set text property to input form
      document.getElementById("text_label").value = originalText
      document.getElementById("text_size").value = originalFontSize

      MicroModal.show("label_edit_modal", {
        onShow: (modal) => {
          document.getElementById("label_edit_ok").onclick = () => {
            targetText.textContent = document.getElementById("text_label").value
            targetText.style.fontSize = document.getElementById("text_size").value
          }
        },
      })
    })
    text.addEventListener("dblclick", (e) => {
      text.style.fill = document.getElementById("colorpicker").style.color
    })
  }

  // Setup SVG pan/zoom
  const panzoom = Panzoom(svg, {
    canvas: true,
    minScale: 0.7,
    maxScale: 10,
  })
  svg.parentElement.addEventListener("wheel", panzoom.zoomWithWheel)
  document.getElementById("zoom_in").addEventListener("click", panzoom.zoomIn)
  document.getElementById("zoom_out").addEventListener("click", panzoom.zoomOut)
  document.getElementById("zoom_reset").addEventListener("click", panzoom.reset)

  // Download PNG Image
  document.getElementById("download_png").addEventListener("click", () => {
    svg.addEventListener(
      "panzoomreset",
      () => {
        downloadPng(svg)
      },
      { once: true }
    )
    panzoom.reset({ animate: false })
  })

  // Download SVG Image
  document.getElementById("download_svg").addEventListener("click", () => {
    svg.addEventListener(
      "panzoomreset",
      () => {
        downloadSvg(svg)
      },
      { once: true }
    )
    panzoom.reset({ animate: false })
  })

  // Setup Features Table
  const featuresTable = new Tabulator("#tabulator", {
    data: Object.values(FEATURES_JSON),
    height: 500,
    layout: "fitData",
    pagination: true,
    paginationSize: 50,
    paginationCounter: "pages",
    paginationSizeSelector: [25, 50, 100, true],
    groupBy: "track",
    groupToggleElement: "header",
    rowContextMenu: [
      {
        label: (row) => {
          return `Search for protein_id='${row.getData().protein_id}' in NCBI Protein database`
        },
        action: (e, row) => {
          const targetUrl = `https://www.ncbi.nlm.nih.gov/protein/${row.getData().protein_id}`
          open(targetUrl, "_blank")
        },
        disabled: (row) => {
          return row.getData().protein_id === "na"
        },
      },
      {
        label: "Search for similar protein sequences using NCBI BLASTP",
        action: (e, row) => {
          const query = encodeURIComponent(convertFeatureRowToFasta(row))
          const baseUrl = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
          const params = `PAGE_TYPE=BlastSearch&PROGRAM=blastp&PAGE=Proteins&QUERY=${query}`
          const targetUrl = `${baseUrl}?${params}`
          open(targetUrl, "_blank")
        },
        disabled: (row) => {
          return row.getData().translation === "na"
        },
      },
      {
        label: "Copy protein fasta to clipboard",
        action: (e, row) => {
          const fasta = convertFeatureRowToFasta(row)
          navigator.clipboard.writeText(fasta).then(
            () => {
              alert(`Copy protein fasta to clipboard.\n\n${fasta}`)
            },
            () => {
              alert("Failed to copy protein fasta to clipboard.")
            }
          )
        },
        disabled: (row) => {
          return row.getData().type !== "CDS" || row.getData().translation === "na"
        },
      },
      {
        label: "Filter by target column field value",
        menu: getFilterMenuList(),
      },
    ],
    columns: [
      { title: "gid", field: "gid", visible: false },
      { title: "translation", field: "translation", visible: false },
      {
        title: "track",
        field: "track",
        headerSort: false,
        minWidth: 150,
        headerFilter: "list",
        headerFilterParams: { valuesLookup: "all" },
      },
      {
        title: "segment",
        field: "segment",
        headerSort: false,
        minWidth: 150,
        headerFilter: "list",
        headerFilterParams: { valuesLookup: "all" },
      },
      {
        title: "location",
        field: "location",
        headerSort: false,
        minWidth: 200,
        headerTooltip: "start - end (strand)",
      },
      {
        title: "length",
        field: "length",
        headerSort: false,
        minWidth: 80,
        formatter: (cell) => {
          return cell.getValue().toLocaleString("en-US")
        },
      },
      {
        title: "type",
        field: "type",
        headerSort: false,
        minWidth: 80,
        headerFilter: "list",
        headerFilterParams: { valuesLookup: "all" },
        headerTooltip: "Feature Type (e.g. CDS, rRNA, tRNA, etc...)",
      },
      {
        title: "gene",
        field: "gene",
        headerSort: false,
        minWidth: 80,
        headerFilter: "list",
        headerFilterParams: { valuesLookup: "all", autocomplete: true },
      },
      {
        title: "pseudo",
        field: "pseudo",
        headerSort: false,
        minWidth: 80,
        headerFilter: "list",
        headerFilterParams: { values: ["true", "false"] },
        headerTooltip: "true/false value of whether feature is a pseudogene or not",
      },
      {
        title: "protein_id",
        field: "protein_id",
        headerSort: false,
        minWidth: 150,
        headerFilter: "list",
        headerFilterParams: { valuesLookup: "all", autocomplete: true },
      },
      {
        title: "product",
        field: "product",
        tooltip: true,
        headerSort: false,
        minWidth: 150,
        headerFilter: "list",
        headerFilterParams: { valuesLookup: "all", autocomplete: true },
      },
    ],
  })
  const gid2FeaturePathList = getGid2FeaturePathList(svgPathList)
  featuresTable.on("rowMouseOver", function (e, row) {
    gid2FeaturePathList[row.getData().gid].forEach((featurePath) => {
      featurePath.insertAdjacentHTML("beforeend", SVG_ANIMATE_HTML)
    })
  })
  featuresTable.on("rowMouseOut", (e, row) => {
    gid2FeaturePathList[row.getData().gid].forEach((featurePath) => {
      featurePath.querySelector("animate").remove()
    })
  })
  featuresTable.on("rowDblClick", (e, row) => {
    gid2FeaturePathList[row.getData().gid].forEach((featurePath) => {
      featurePath.style.fill = document.getElementById("colorpicker").style.color
    })
  })
  featuresTable.on("menuOpened", (row) => {
    row.select()
  })
  featuresTable.on("menuClosed", (row) => {
    row.deselect()
  })
  document.getElementById("apply_color").addEventListener("click", () => {
    const activeTableRows = featuresTable.getRows("active")
    activeTableRows.forEach((row) => {
      gid2FeaturePathList[row.getData().gid].forEach((featurePath) => {
        featurePath.style.fill = document.getElementById("colorpicker").style.color
      })
    })
  })
  document.getElementById("clear_filter").addEventListener("click", () => {
    featuresTable.clearHeaderFilter()
  })
  document.getElementById("download_table").addEventListener("click", () => {
    featuresTable.download("csv", "features_table.tsv", { delimiter: "\t" })
  })
})

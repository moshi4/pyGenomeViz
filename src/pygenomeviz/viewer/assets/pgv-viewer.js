/**
 * Convert feature ID to proper label for tooltip display
 * @param {string} featureId
 * @returns {string} Feature info for tooltip display
 */
function convertFeatureIdToLabel(featureId) {
  const matchResult = featureId.match(/Feature_(\d+)_(\d+)_(-*\d+)_type_(.+)_gene_(.+)_protein_id_(.+)_product_(.+)/)
  let [start, end, strand, type, gene, proteinId, product] = matchResult.slice(1, 8)
  const locale = "en-US"
  length = Number(end - start).toLocaleString(locale)
  start = Number(start).toLocaleString(locale)
  end = Number(end).toLocaleString(locale)
  strand = strand === "-1" ? "-" : "+"
  product = product.replaceAll("_", " ")
  let label = `location: ${start} - ${end} (${strand})\nlength: ${length} bp`
  if (type !== "na") {
    label += `\ntype: ${type}\ngene: ${gene}\nprotein_id: ${proteinId}\nproduct: ${product}`
  }
  return label
}
/**
 * Convert link ID to proper label for tooltip display
 * @param {string} linkId
 * @returns {string} Link info for tooltip display
 */
function convertLinkIdToLabel(linkId) {
  const matchResult = linkId.match(/Link_(.+)_(\d+)_(\d+)_(.+)_(\d+)_(\d+)_(.+)$/)
  let [n1, s1, e1, n2, s2, e2, ident] = matchResult.slice(1, 8)
  const locale = "en-US"
  s1 = Number(s1).toLocaleString(locale)
  e1 = Number(e1).toLocaleString(locale)
  s2 = Number(s2).toLocaleString(locale)
  e2 = Number(e2).toLocaleString(locale)
  let label = `1. ${n1} (${s1} - ${e1} bp)\n2. ${n2} (${s2} - ${e2} bp)`
  if (ident !== "na") {
    label += `\nIdentity: ${ident}%`
  }
  return label
}
/**
 * Save as PNG image
 * @param {HTMLElement} svgNode
 * @param {string} fileName
 */
function saveAsPng(svgNode, fileName = "image.png") {
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
  image.src = "data:image/svg+xml;charset=utf-8;base64," + btoa(decodeURIComponent(encodeURIComponent(svgData)))
}
/**
 * Save as  SVG image
 * @param {HTMLElement} svgNode
 * @param {string} fileName
 */
function saveAsSvg(svgNode, fileName = "image.svg") {
  const svgData = new XMLSerializer().serializeToString(svgNode)
  const svgBlob = new Blob([svgData], { type: "image/svg+xml;charset=utf-8" })

  const a = document.createElement("a")
  a.href = URL.createObjectURL(svgBlob)
  a.setAttribute("download", fileName)
  a.dispatchEvent(new MouseEvent("click"))
}

$(document).ready(function () {
  // Remove unnecessary previous saved HTML elements
  document.querySelectorAll(".ui-tooltip").forEach((e) => e.remove())
  document.querySelectorAll(".sp-container").forEach((e) => e.remove())
  document.querySelectorAll(".ui-helper-hidden-accessible").forEach((e) => e.remove())

  const svg = document.querySelector("#svg_canvas>svg")

  // Set colorpicker
  const initColor = "red"
  document.getElementById("colorpicker").style.color = initColor
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
      if (pickColor === "#ffffff") {
        pickColor = "transparent"
      }
      document.getElementById("colorpicker").style.color = pickColor
    },
  })

  const allPaths = svg.getElementsByTagName("path")
  for (let path of allPaths) {
    // Set feature & link tooltip
    const pathId = path.parentNode.id
    if (!pathId.startsWith("Feature") && !pathId.startsWith("Link")) {
      continue
    }

    let $pathObj = $("#" + pathId + ">path")
    $pathObj.attr("title", "")
    let tooltipLabel = ""
    if (pathId.startsWith("Feature")) {
      tooltipLabel = convertFeatureIdToLabel(pathId)
    } else if (pathId.startsWith("Link")) {
      tooltipLabel = convertLinkIdToLabel(pathId)
    }
    $pathObj.tooltip({ content: tooltipLabel, show: false, hide: false, track: true })

    // Set picked color as facecolor
    path.addEventListener("dblclick", () => {
      path.style.fill = document.getElementById("colorpicker").style.color
    })

    // Highlight selected path objects
    const originalStroke = path.style.stroke
    const originalStrokeWidth = path.style["stroke-width"]
    path.addEventListener("mouseover", () => {
      path.style.stroke = "black"
      path.style["stroke-width"] = "1.0"
    })
    path.addEventListener("mouseout", () => {
      path.style.stroke = originalStroke
      path.style["stroke-width"] = originalStrokeWidth
    })
  }

  // Edit Label Text
  const allTexts = svg.getElementsByTagName("text")
  for (let text of allTexts) {
    // Show edit label text & size dialog
    text.addEventListener("contextmenu", (e) => {
      e.preventDefault()
      const originalText = text.textContent
      const originalFont = parseInt(text.style.font)
      $("#text_label").val(originalText)
      $("#text_size").val(originalFont)
      $("#label_edit_dialog").dialog({
        modal: true,
        height: 250,
        width: 400,
        title: "Edit Label Text",
        dialogClass: "no-close",
        buttons: {
          Cancel: function () {
            $(this).dialog("close")
          },
          OK: function () {
            text.textContent = $("#text_label").val()
            text.style.font = $("#text_size").val() + "px 'sans-serif'"
            $(this).dialog("close")
          },
        },
      })
    })
    // Set picker color as label color
    text.addEventListener("dblclick", () => {
      text.style.fill = document.getElementById("colorpicker").style.color
    })
  }

  // Help Manual Dialog
  const helpIcon = document.getElementById("help_icon")
  helpIcon.addEventListener("click", () => {
    $("#help_dialog").dialog({
      modal: false,
      height: 600,
      width: 700,
      title: "Help Manual",
      dialogClass: "no-close",
      buttons: {
        Close: function () {
          $(this).dialog("close")
        },
      },
    })
  })

  // SVG Pan & Zoom setting
  const panzoom = Panzoom(svg, {
    canvas: true,
    minScale: 0.7,
    maxScale: 5,
  })
  svg.parentElement.addEventListener("wheel", panzoom.zoomWithWheel)
  document.getElementById("zoom_in").addEventListener("click", panzoom.zoomIn)
  document.getElementById("zoom_out").addEventListener("click", panzoom.zoomOut)
  document.getElementById("zoom_reset").addEventListener("click", panzoom.reset)

  // Save as PNG Image
  document.getElementById("png_save").addEventListener("click", () => {
    svg.addEventListener(
      "panzoomreset",
      () => {
        saveAsPng(svg)
      },
      { once: true }
    )
    panzoom.reset({ animate: false })
  })

  // Save as SVG Image
  document.getElementById("svg_save").addEventListener("click", () => {
    svg.addEventListener(
      "panzoomreset",
      () => {
        saveAsSvg(svg)
      },
      { once: true }
    )
    panzoom.reset({ animate: false })
  })
})

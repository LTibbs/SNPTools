//Define gloabl variables
let main_id = "";
let main_option = "both";
let can_flag = false;
let gene_model, transcript, protein, gm_can, protein_can, gs, gn, unip_id, unip_desc, gene_model_file, pfamArray, alignment_length;
let gene_model_length = 0;
let pfam_domain_array = [];
let view_type = "Species";

let X2_to_X_array = {}; //converts reference coordinates to MSA coordinates
let GM_array = {};
let GN_array = {};
let X_to_WT = {}; //converts  MSA coordinates to wild-typoe amino acid
let X_to_X2 = {}; //converts  MSA coordinates to reference coordinates
let GN_size = 0; //Gene number

let scaleFactor = 0;
let scaleFactorPan = 0;
const window_length = 1200;  //View window size, everything is normalized based on this value
let scaleFactorZoom = window_length / 50;  //Normalization of the size of each cell in the heatmap
let scaleFactorZoomPan = window_length / 50;  //Normalization of the size of each cell in the heatmap
let currentlyVisibleTooltip = null;  //Controls what tooltip is visible

var wgs_status = true;
var wgs2024_status = false;
var wgs2026_status = false;

//This section loads the MaizeGDB specific annotations including Uniprot and gene annotations, this section would be customized for other model specific annotations
let uniprot_filename = './uniprot/' + gene_model + '.tsv';
let synonym_filename = './synonym/maize_synonym.tsv';

/* ---------------- ESM model toggle ---------------- */

// Directory names exactly as you described:
const ESM_MODELS = ["ESM1", "ESM2", "ESM3"];
let currentESM = "ESM2"; // default

function getESMFromUrlOrStorage() {
  // support.js already provides getQueryStringValue() (you use it in index.html) :contentReference[oaicite:6]{index=6}
  const fromUrl = (typeof getQueryStringValue === "function") ? getQueryStringValue("esm") : null;
  const fromStorage = localStorage.getItem("panEffectESM");

  const candidate = (fromUrl || fromStorage || "ESM2").toUpperCase();
  return ESM_MODELS.includes(candidate) ? candidate : "ESM2";
}

function setESMModel(newModel) {
  const m = (newModel || "").toUpperCase();
  currentESM = ESM_MODELS.includes(m) ? m : "ESM2";
  localStorage.setItem("panEffectESM", currentESM);

  // keep links shareable: add/update ?esm=ESM#
  try {
    const u = new URL(window.location.href);
    u.searchParams.set("esm", currentESM);
    history.replaceState({}, "", u.toString());
  } catch (e) {
    // ignore if URL API not available
  }
}

function csvFileForProtein(proteinId) {
  return `./csv/${currentESM}/${proteinId}.csv`;
}

function heatmapFileForCanonicalTranscript(gmCan) {
  return `./heatmap/${currentESM}/${gmCan}.tsv`;
}

/* optional: remember last flags so a toggle can re-render */
let last_b73_flag = true;
let last_pan_flag = true;

function clearHeatmapViews() {
  // Clear containers so re-render doesn’t stack DOM
  const idsToClear = [
    "full-heatmap",
    "zoomed-heatmap",
    "heatNumberLine",
    "zoomNumberLine",
    "zoomWTLine",
    "heatNumberLine-pan",
    "zoomNumberLine-pan",
    "zoomed-heatmap-pan"
  ];
  idsToClear.forEach(id => {
    const el = document.getElementById(id);
    if (el) el.innerHTML = "";
  });

  // highlight boxes are re-added by your render code; keep simple:
  const full = document.getElementById("full-heatmap");
  if (full) {
    const hb = document.createElement("div");
    hb.id = "highlight-box";
    full.appendChild(hb);
  }
}

function initESMSelector() {
  currentESM = getESMFromUrlOrStorage();
  setESMModel(currentESM);

  const sel = document.getElementById("esmModelSelect");
  if (!sel) return;

  sel.value = currentESM;

  sel.addEventListener("change", async (e) => {
    onLoadingIcon();
    setESMModel(e.target.value);

    // Re-load score-driven views from the new directories
    clearHeatmapViews();

    // IMPORTANT: openGM sets the global gene_model_file used for the download link :contentReference[oaicite:7]{index=7}
    const ok = await openGM(protein);
    if (ok) {
      await fetchDataAndSetup(protein, protein_can, last_b73_flag, last_pan_flag);
    }
    offLoadingIcon();
  });
}


async function openUniprot(id) {
    try {
        let uniprot_filename = './uniprot/' + id + '.tsv';

        let response = await fetch(uniprot_filename);
        if (!response.ok) {
            return false;
        }

        let data = await response.text();

        if(data.includes("<html") || data.includes("<script"))
        {
            return false;
        }

        let rows = data.split('\n');

        rows.forEach(row => {
            let dataU = row.split("\t");
            if (dataU.length >= 7) {
                let [gm, gm_can_local, unip_id_local, uni_desc_local, goU, gs_local, gn_local] = dataU;
                gm_can = gm_can_local;
                protein_can = gm_can.replace("_T", "_P");
                gs = gs_local;
                gn = gn_local;
                unip_id = unip_id_local;
                unip_desc = uni_desc_local;
                if (can_flag) {
                    transcript = gm_can;
                    protein = gm_can.replace("_T", "_P");
                }
            }
        });

        return true;
    } catch (error) {
        console.error('There was a problem with the fetch operation:', error.message);
        return false;
    }
}

//This function loads the CSV containing the positions of the heatmap, the variant effect score, and the amino acids for the WT and Substitution
//HEADER in CSV: X,Y,Score,WT,Sub
//SMAPLE DATA: 1,9,-9.54,M,K

async function openGM(protein) {
gene_model_file = csvFileForProtein(protein);
gene_model_length = 0;

try {
    let response = await fetch(gene_model_file);

    if (!response.ok) {
        console.error("Failed to fetch gene model CSV file.");
        return false;
    }

    let data = await response.text();

    if(data.includes("<html") || data.includes("<script"))
    {
        return false;
    }

    let rows = data.trim().split("\n");
    let lastLine = rows[rows.length - 1];

    if (lastLine) {
        let lastData = lastLine.split(",");

        if (lastData.length > 0) {
            gene_model_length = parseInt(lastData[0], 10);
            scaleFactor = window_length / gene_model_length;
        }
    }

    return true;

} catch (error) {
    console.error("There was an error:", error.message);
    return false;  // Returns -1 to indicate an error occurred, but this can be adjusted to your needs
}
}

async function fetchAndProcessQueryData(queryFilename) {
    let queryArray = [];

    try {
        // Fetch the TSV content
        let response = await fetch(queryFilename);
        if (!response.ok) {
            document.getElementById("pan-genome").innerHTML = '<div class="text-content" style="color: red;">Pan-genome data for the gene, transcript, or protein id <b>' + main_id + '</b> was not found.  Use the search tab to enter a new search term.</div>';
            return false;
            //throw new Error("Failed to fetch the file.");
        }

        let data = await response.text();

        if(data.includes("<html") || data.includes("<script"))
        {
            document.getElementById("pan-genome").innerHTML = '<div class="text-content" style="color: red;">Pan-genome data for the gene, transcript, or protein id <b>' + main_id + '</b> was not found.  Use the search tab to enter a new search term.</div>';
            return false;
        }

        const lines = data.trim().split("\n");

        // Loop through the lines (skipping the header) to process the data
        for (let i = 1; i < lines.length; i++) {
            const [X, X2, WT] = lines[i].split("\t");

            // Populating the X_to_X2 and X_to_WT maps
            X_to_X2[X] = X2;
            X_to_WT[X] = WT;

            // Populating the queryArray and other variables
            queryArray.push({
                X: parseInt(X, 10),
                X2: X2,
                WT: WT
            });

            if (X !== null && parseInt(X, 10) > 0) {
                alignment_length = parseInt(X, 10);
                scaleFactorPan = window_length / alignment_length;
            }
            X2_to_X_array[X2] = parseInt(X, 10);
        }

        return true;

    } catch (error) {
        document.getElementById("pan-genome").innerHTML = '<div class="text-content" style="color: red;">Pan-genome data for the gene, transcript, or protein id <b>' + main_id + '</b> was not found.  Use the search tab to enter a new search term.</div>';
        console.log("Error fetching or processing the file:", error);
        return false;
    }
}

async function fetchTargetData(targetFilename) {
    try {
        // Fetch the TSV content
        let response = await fetch(targetFilename);
        let data = await response.text();

        // Split the TSV content into lines
        const lines = data.trim().split("\n");

        // Loop through the lines (skipping the header) to process the data
        for(let i = 1; i < lines.length; i++) {
            const [Y, GM, GN] = lines[i].split("\t");
            GM_array[Y] = GM;
            GN_array[Y] = GN;
        }

        GN_size = lines.length;
        let heatmapContainer = document.querySelector('.heatmap-container-pan');
        let heatmap_div_height = (GN_size  * 5);
        heatmapContainer.style.height = `${heatmap_div_height}px`;

        const zoomedContainer = document.getElementById('heatmap-container-zoom-pan');
        //zoomedContainer.style.height = ((GN_size - 1) * 20) + "px";
        zoomedContainer.style.height = ((GN_size) * 20) + "px";
        return true;

    } catch (error) {
        document.getElementById("pan-genome").innerHTML = '<div class="text-content" style="color: red;">Pan-genome data for the gene, transcript, or protein id <b>' + main_id + '</b> was not found.  Use the search tab to enter a new search term.</div>';
        console.log("Target Error fetching or processing the TSV file:", error);
    }
}

async function fetchPfamData(transcript) {
    let pfam_filename = './pfam/B73/' + transcript + '.tsv';
    let query_filename = './query/' + gm_can + '.tsv';

    // Reset PFAM state so it can’t accumulate across reloads
    pfam_domain_array = [];
    pfamArray = [];

    try {
        await fetchAndProcessQueryData(query_filename);
        populateSummary();
        loadAndDisplayTraits(gene_model);

        let response = await fetch(pfam_filename);

        if (!response.ok) {
            offLoadingIcon();
            document.getElementById("pfam-wrap").innerHTML = '<div class="text-content"><br>There are no PFam domains for this gene model.</div>';
            document.getElementById("pfam-wrapl-pan").innerHTML = '<div class="text-content"><br>There are no PFam domains for this gene model.</div>';
            return;
        }

        let data = await response.text();
        let rows = data.trim().split("\n");

        rows.forEach(row => {
            let columns = row.split("\t");

            // Assuming each line has 13 or 14 columns
            if (columns.length >= 13) {
                let gm, gm_length, pfam_name, pfam_desc, pfam_start, pfam_end, pfam_score, inter_name, inter_desc, go_terms;

            if (columns.length == 13)
            {
                [gm, , gm_length, , pfam_name, pfam_desc, pfam_start, pfam_end, pfam_score, , , inter_name, inter_desc] = columns;
                go_terms = "N/A";
            }

            if (columns.length == 14)
            {
                [gm, , gm_length, , pfam_name, pfam_desc, pfam_start, pfam_end, pfam_score, , , inter_name, inter_desc, go_terms] = columns;
            }
                pfam_name = pfam_name || "N/A";
                pfam_desc = pfam_desc || "N/A";
                inter_name = inter_name || "N/A";
                inter_desc = inter_desc || "N/A";
                go_terms = go_terms || "N/A";

                let domainObject = {
                    start: parseInt(pfam_start, 10),
                    end: parseInt(pfam_end, 10),
                    start2: parseInt(X2_to_X_array[pfam_start], 10),
                    end2: parseInt(X2_to_X_array[pfam_end], 10),
                    pfam_id: pfam_name,
                    pfam_name: pfam_desc,
                    inter_id: inter_name,
                    inter_name: inter_desc,
                    go: go_terms
                };

                pfam_domain_array.push(domainObject);
            }
        });

        offLoadingIcon();
        return pfam_domain_array;

    } catch (error) {
        offLoadingIcon();
        console.log("Error in fetchPfamData:", error);
    }
}

// ===== Native SNPTools entry point =====================================
// The engine sets the shared globals (main_id, main_option, gene_model,
// transcript, protein, can_flag, currentESM) then calls runPanEffect().
// No URL, tab, or localStorage coupling — the panel drives view selection.
async function runPanEffect() {

    onLoadingIcon();

    try {
        const result_syn = await checkSynonym(gene_model);
        if (result_syn !== false) { gene_model = result_syn; }
    } catch (error) {
        console.log("No synonym found");
    }

    try {
        const UniProtresult = await openUniprot(gene_model);
        if (!UniProtresult) {
            if (main_id) { errorInnerHTML(); } else { emptyInnerHTML(); }
            offLoadingIcon();
            return false;
        }

        const result = await openGM(protein);
        scaleFactor = window_length / gene_model_length;

        // main_option is 'both' | 'b73' | 'pan'; pan-genome is canonical-only
        let b73_flag = (main_option !== "pan");
        let pan_flag = (main_option !== "b73");
        if (protein !== protein_can) { pan_flag = false; }

        last_b73_flag = b73_flag;
        last_pan_flag = pan_flag;

        if (result) {
            await fetchDataAndSetup(protein, protein_can, b73_flag, pan_flag);
            return true;
        } else {
            if (main_id) { errorInnerHTML(); }
            offLoadingIcon();
            return false;
        }
    } catch (error) {
        console.error("[PanEffect] run error:", error);
        offLoadingIcon();
        return false;
    }
}

//This function checks the synonym file for gene names
async function checkSynonym(id) {
  return new Promise((resolve, reject) => {
    // Use the Fetch API to fetch the file content
    fetch(synonym_filename)
      .then((response) => {
        if (!response.ok) {
          throw new Error(`Failed to fetch file: ${response.statusText}`);
        }
        return response.text(); // Get the file content as text
      })
      .then((fileContents) => {
        // Split the file contents into lines
        const lines = fileContents.split('\n');
        // Iterate through each line to check for the ID in the 2nd column
        let xx = 1;
        for (const line of lines) {
          const columns = line.split(/\s+|\t+/);
          if (columns.length >= 2 && columns[1] == id) {
            //console.log(`ID '${id}' found.`);
            resolve(columns[0]); // Resolve the Promise with the value
            return; // Exit the function after finding the ID
          }
        }
        // If the ID is not found, reject the Promise with false
        reject(false);
      })
      .catch((error) => {
        reject(error);
      });
  });
}

//Main function that loads most of the data sets
async function fetchDataAndSetup(protein, protein_can, b73_flag, pan_flag) {
    try {
        let gene_model_file = csvFileForProtein(protein);
        let target_filename = './target/' + gm_can + '.tsv';
        let query_filename = './query/' + gm_can + '.tsv';
        let pfam_filename = './pfam/B73/' + transcript + '.tsv';
        let pfam_flag = true;
        let response = await fetch(gene_model_file);

        if (!response.ok) {
            throw new Error("Failed to fetch gene model CSV file.");
        }

        let responseP = await fetch(pfam_filename);
        let dataP = await responseP.text();

        if (!responseP.ok) {
            pfam_flag = false;
        }

        if(dataP.includes("<html") || dataP.includes("<script"))
        {
            pfam_flag = false;
        }

        let data = await response.text();

        const parsedData = d3.csvParse(data);
        if(b73_flag)
        {
            updateHeatmapZoom(parsedData, 1);
            updateHighlightBox(1); // was `+this.value`, which is meaningless here (not inside a slider event handler) and evaluated to NaN
            createZoomedWTLine(document.getElementById('zoomWTLine'), parsedData, 1);
            renderHeatmap(parsedData);
        }

        if(pan_flag)
        {
            let responseQ = await fetchAndProcessQueryData(query_filename);

            if (!responseQ) {
                pan_flag = false;
            }
        }

        let responseTarget;
        let parsedDataCan;

        if(b73_flag)
        {
            sliderInnerHTML("slider");

            let heatmapContainer2 = document.querySelector('.heatmap-container');
            let heatmap_div_height2 = (30  * 10);
            heatmapContainer2.style.height = heatmap_div_height2 + 'px';

            const zoomedContainer2 = document.getElementById('heatmap-container-zoom');
            zoomedContainer2.style.height = ((29) * 14) + "px";
        }

        if(pan_flag)
        {
            try {
                    responseTarget = await fetchTargetData(target_filename);
                    if (responseTarget) {

                        sliderInnerHTMLPan("slider-pan");

                        let heatmap_filename = heatmapFileForCanonicalTranscript(gm_can);

                        let response_can = await fetch(heatmap_filename);

                        if (!response_can.ok) {
                            document.getElementById("pan-genome").innerHTML = '<div class="text-content" style="color: red;">Pan-genome data for the gene, transcript, or protein id <b>' + main_id + '</b> was not found.  Use the search tab to enter a new search term.</div>';
                            throw new Error("Failed to fetch gene model canonical CSV file.");
                        }
                        let data_can = await response_can.text();

                        if(data_can.includes("<html") || data_can.includes("<script"))
                        {
                            document.getElementById("pan-genome").innerHTML = '<div class="text-content" style="color: red;">Pan-genome data for the gene, transcript, or protein id <b>' + main_id + '</b> was not found.  Use the search tab to enter a new search term.</div>';
                            throw new Error("Failed to fetch gene model canonical CSV file.");
                        }

                        parsedDataCan = d3.tsvParse(data_can);
    

                        updateHeatmapZoomPan(parsedDataCan, 1);
                        updateHighlightBoxPan(1);
                        renderHeatmapPan(parsedDataCan);
                    } else {
                        pan_flag = false;
                    }
                } catch (error) {
                    offLoadingIcon();
                    document.getElementById("pan-genome").innerHTML = '<div class="text-content" style="color: red;">Pan-genome data for the gene, transcript, or protein id <b>' + main_id + '</b> was not found.  Use the search tab to enter a new search term.</div>';
                    console.error("There was an error:", error.message);
                }
      }

        fetchPfamData(transcript).then(pfam_domain_array => {
            pfamArray = pfam_domain_array;
            scaleFactor = window_length / gene_model_length;

            if(b73_flag)
            {
                if (pfam_flag) {
                    renderDomains(document.getElementById('pfamGeneModel'), pfamArray);
                    createNumberLine(document.getElementById('pfamNumberLine'));
                } else {
                    document.getElementById("pfam-wrap").innerHTML = '<div class="text-content"><br>There are no PFam domains for this gene model.</div>';
                }

                createNumberLine(document.getElementById('heatNumberLine'));
                // Read the slider's CURRENT value rather than hardcoding 1: this callback
                // is not awaited, so it can resolve before OR after centerOnVariant() moves
                // the window to an incoming variant's position. Hardcoding 1 here silently
                // reverted the zoomed-region ticks back to the default range whenever this
                // fetch happened to finish second, even though the heatmap cells/WT line
                // underneath were already showing the right region.
                {
                    const zoomSliderEl = document.getElementById('zoom-slider');
                    const curStart = zoomSliderEl ? (+zoomSliderEl.value || 1) : 1;
                    createZoomedNumberLine(document.getElementById('zoomNumberLine'), curStart);
                }
                loadDSSP();
            }

            if(pan_flag)
            {
                if (pfam_flag) {
                    renderDomainsPan(document.getElementById('pfamGeneModel-pan'), pfamArray);
                    createNumberLinePan(document.getElementById('pfamNumberLine-pan'));
                } else {
                    document.getElementById("pfam-wrap-pan").innerHTML = '<div class="text-content"><br>There are no PFam domains for this gene model.</div>';
                }

                createNumberLinePan(document.getElementById('heatNumberLine-pan'));
                // same fix as the B73 side above, for the pan-genome slider
                {
                    const zoomSliderPanEl = document.getElementById('zoom-slider-pan');
                    const curStartPan = zoomSliderPanEl ? (+zoomSliderPanEl.value || 1) : 1;
                    createZoomedNumberLinePan(document.getElementById('zoomNumberLine-pan'), curStartPan);
                }
                loadDSSPPan();
            }
        });

        if(b73_flag)
        {
            // Add listener for slider change
            const slider = document.getElementById("zoom-slider");
            const sliderValueDisplay = document.getElementById("slider-value");

            // Function to handle the change event
            function handleRadioChange(event) {
                wgs_status = (event.target.value === 'all');
                wgs2024_status = (event.target.value === 'maize2024');
                wgs2026_status = (event.target.value === 'maize2026');

                const slider = document.getElementById("zoom-slider");
                const sliderValueDisplay = document.getElementById("slider-value");

                // Update heatmap based on slider value
                renderHeatmap(parsedData);
                updateHeatmapZoom(parsedData, +slider.value);

            }

            // Add event listeners to the radio buttons
            document.querySelectorAll('input[name="variantEffect"]').forEach((radio) => {
                radio.addEventListener('change', handleRadioChange);
                //console.log("Change");
            });

            slider.addEventListener("input", function() {
                // Update the display value
                sliderValueDisplay.textContent = this.value;

                // Update heatmap based on slider value
                updateHeatmapZoom(parsedData, +this.value);
                createZoomedWTLine(document.getElementById('zoomWTLine'), parsedData, +this.value);
                updateHighlightBox(+this.value);
                createZoomedNumberLine(document.getElementById('zoomNumberLine'), +this.value);
            });
        }

        if(pan_flag)
        {
            // Add listener for slider change
            const sliderPan = document.getElementById("zoom-slider-pan");
            const sliderValueDisplayPan = document.getElementById("slider-value-pan");

            sliderPan.addEventListener("input", function() {
                // Update the display value
                sliderValueDisplayPan.textContent = this.value;

                // Update heatmap based on slider value
                updateHeatmapZoomPan(parsedDataCan, +this.value);
                updateHighlightBoxPan(+this.value);
                createZoomedNumberLinePan(document.getElementById('zoomNumberLine-pan'), +this.value);
            });
        }

    } catch (error) {
        offLoadingIcon();
        console.error("There was an error:", error.message);
    }
}

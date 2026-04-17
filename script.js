let stage, component;
let spin=false;
let ligandRep=null;
let activeRep=null;

// Protein name → PDB
const proteinMap = {
"hemoglobin":"6HHB",
"insulin":"1INS",
"lysozyme":"1LYZ",
"myoglobin":"1MBN",
"dna":"1BNA",
"spike":"6VSB",
"covid":"6LU7",
"cas9":"5XNL",
"ribosome":"3J3Q"
};

// INIT
function init(){

stage = new NGL.Stage("viewer");

// click residue info
stage.signals.clicked.add(function(p){

if(p && p.atom){
let a = p.atom;

residueInfo.innerHTML = `
<b>Residue Info</b><br>
Name: ${a.resname}<br>
Chain: ${a.chainname}<br>
Position: ${a.resno}
`;
}

});
}

// LOAD PROTEIN
function loadProtein(){

document.getElementById("viewer").innerHTML="";
init();

let input = searchProtein.value.toLowerCase().trim();
let pdb = proteinMap[input] || input || "1CRN";

stage.loadFile("rcsb://" + pdb).then(comp=>{

component = comp;

// base structure
component.addRepresentation(rep.value,{colorScheme:"chainname"});
component.addRepresentation("cartoon",{colorScheme:"sstruc"});

component.autoView();

fetchInfo(pdb);

}).catch(()=>{
alert("Protein not found");
});
}

// ROTATE
function toggleSpin(){
spin=!spin;
stage.setSpin(spin);
}

// RESET
function resetView(){
component.autoView();
stage.setSpin(false);
}

// LIGAND
function toggleLigand(){

if(!component) return;

if(!ligandRep){
ligandRep = component.addRepresentation("ball+stick",{
sele:"hetero and not water",
color:"red"
});
}else{
component.removeRepresentation(ligandRep);
ligandRep=null;
}
}

// ACTIVE SITE
function toggleActiveSite(){

if(!component) return;

if(!activeRep){
activeRep = component.addRepresentation("spacefill",{
sele:"within 4 of (hetero and not water)",
color:"yellow",
opacity:0.7
});
}else{
component.removeRepresentation(activeRep);
activeRep=null;
}
}

// ANALYSIS
function analyzeProtein(){

if(!component){
alert("Load protein first");
return;
}

let helix=0, sheet=0, coil=0, total=0;
let residues=0, atoms=0, ligands=0;
let amino={}, chains={};

component.structure.eachResidue(r=>{

total++; residues++;

chains[r.chainname]=(chains[r.chainname]||0)+1;
amino[r.resname]=(amino[r.resname]||0)+1;

if(r.isHet) ligands++;

if(r.sstruc==="h") helix++;
else if(r.sstruc==="e") sheet++;
else coil++;
});

component.structure.eachAtom(()=>atoms++);

let h=(helix/total*100).toFixed(1);
let s=(sheet/total*100).toFixed(1);
let c=(coil/total*100).toFixed(1);

let chainInfo = Object.entries(chains)
.map(x=>`Chain ${x[0]}: ${x[1]}`)
.join("<br>");

let topAA = Object.entries(amino)
.sort((a,b)=>b[1]-a[1])
.slice(0,5)
.map(x=>`${x[0]} (${x[1]})`)
.join("<br>");

analysisPanel.innerHTML=`
<h3>Protein Analysis</h3>

<b>Secondary Structure</b><br>
Helix: ${h}%<br>
Sheet: ${s}%<br>
Coil: ${c}%<br>

<br><b>Stats</b><br>
Residues: ${residues}<br>
Atoms: ${atoms}<br>
Ligands: ${ligands}

<br><br><b>Chains</b><br>
${chainInfo}

<br><br><b>Top Amino Acids</b><br>
${topAA}
`;
}

// FETCH INFO
async function fetchInfo(pdb){

try{
let res = await fetch(`https://data.rcsb.org/rest/v1/core/entry/${pdb}`);
let data = await res.json();

proteinInfo.innerHTML=`
<b>${pdb}</b><br>
${data.struct.title}
`;

}catch{
proteinInfo.innerText="No info available";
}
}

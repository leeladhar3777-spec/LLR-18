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

// click residue
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

component.addRepresentation(rep.value,{colorScheme:"chainname"});
component.addRepresentation("cartoon",{colorScheme:"sstruc"});

component.autoView();

fetchInfo(pdb);
fetchSimilarProteins(pdb);

}).catch(()=>{
alert("Protein not found");
});
}

// CONTROLS
function toggleSpin(){ spin=!spin; stage.setSpin(spin); }
function resetView(){ component.autoView(); stage.setSpin(false); }

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

let helix=0,sheet=0,coil=0,total=0;
let residues=0,atoms=0,ligands=0;
let amino={},chains={};

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
.map(x=>`Chain ${x[0]}: ${x[1]}`).join("<br>");

let topAA = Object.entries(amino)
.sort((a,b)=>b[1]-a[1])
.slice(0,5)
.map(x=>`${x[0]} (${x[1]})`).join("<br>");

analysisPanel.innerHTML=`
<h3>Protein Analysis</h3>

Helix: ${h}%<br>
Sheet: ${s}%<br>
Coil: ${c}%<br>

<br>Residues: ${residues}
<br>Atoms: ${atoms}
<br>Ligands: ${ligands}

<br><br><b>Chains</b><br>${chainInfo}

<br><br><b>Top Amino Acids</b><br>${topAA}
`;
}

// INFO
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

// SIMILAR PROTEINS
function fetchSimilarProteins(pdb){

let list=["1CRN","1LYZ","1MBN","6LU7","6VSB"];

let html="<h3>🧬 Similar Proteins</h3>";

list.forEach(id=>{
html += `<button onclick="loadSimilar('${id}')">${id}</button>`;
});

similarProteins.innerHTML = html;
}

function loadSimilar(pdb){
searchProtein.value = pdb;
loadProtein();
}

// ================= ALIGNMENT =================
function alignSequences(){

let s1 = seq1.value.trim().toUpperCase();
let s2 = seq2.value.trim().toUpperCase();

if(!s1 || !s2){
alert("Enter sequences");
return;
}

const match=2, mismatch=-1, gap=-2;

let m=s1.length, n=s2.length;
let dp=Array(m+1).fill().map(()=>Array(n+1).fill(0));

for(let i=0;i<=m;i++) dp[i][0]=i*gap;
for(let j=0;j<=n;j++) dp[0][j]=j*gap;

for(let i=1;i<=m;i++){
for(let j=1;j<=n;j++){
let d = dp[i-1][j-1] + (s1[i-1]===s2[j-1]?match:mismatch);
let u = dp[i-1][j] + gap;
let l = dp[i][j-1] + gap;
dp[i][j]=Math.max(d,u,l);
}
}

// traceback
let a1="", a2="";
let i=m,j=n;

while(i>0 && j>0){
if(dp[i][j] === dp[i-1][j-1] + (s1[i-1]===s2[j-1]?match:mismatch)){
a1=s1[i-1]+a1; a2=s2[j-1]+a2; i--; j--;
}
else if(dp[i][j] === dp[i-1][j] + gap){
a1=s1[i-1]+a1; a2="-"+a2; i--;
}else{
a1="-"+a1; a2=s2[j-1]+a2; j--;
}
}

while(i>0){ a1=s1[i-1]+a1; a2="-"+a2; i--; }
while(j>0){ a1="-"+a1; a2=s2[j-1]+a2; j--; }

let matchLine="";
for(let k=0;k<a1.length;k++){
matchLine += (a1[k]===a2[k]) ? "|" : " ";
}

alignmentResult.innerHTML=`
Score: ${dp[m][n]}

${a1}
${matchLine}
${a2}
`;
    }

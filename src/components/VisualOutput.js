import MoleculeStructure from "./MoleculeStructure";

const VisualOutput = ({ VisualData, SMILES_LIST, subSMILES_LIST }) => {
    return (
      <div>
        {VisualData ? (
          <div className="VisOutputText" style={{
            fontFamily: "monospace",
            position: "absolute",
            top: "60vh",
            left: "30vw",
            fontSize: "30px",
            color: "black"
          }}>
            Generated Molecules
          </div>
        ) : null}
  
        {SMILES_LIST && SMILES_LIST.length > 0 ? (
          <div style={{
            top: "65vh",
            left: "30vw",
            position: "absolute",
            maxHeight: "400px",
            width: "100%"
          }}>
            <div
            id="structure-list"
            className="columns is-desktop"
            style={{ margin: "15px", overflowX: "auto", whiteSpace: "nowrap" }}
            >
              {SMILES_LIST.map((smiles, index) => (
                <div className="column" key={smiles + index}>
                  <MoleculeStructure
                  id={smiles}
                  structure={smiles}
                  height={250}
                  width={250}
                  subStructure={subSMILES_LIST[index]}
                  svgMode
                  />
                </div>
              ))}
            </div>
          </div>
        ) : null}
      </div>
    );
  };
  
  export default VisualOutput;
  
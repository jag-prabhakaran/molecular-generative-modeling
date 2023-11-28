import MoleculeStructure from "./MoleculeStructure";

const VisualOutput = ({ VisualData, SMILES_LIST, subSMILES_LIST }) => {
  const moleculeWidth = 300; // Set a default width for molecules
  const containerWidth = 90; // Set a default container width percentage

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
          width: `${containerWidth}%`,
          maxWidth: "1200px",
          margin: "0 auto",
          overflowX: "auto",
        }}>
          <div
            id="structure-list"
            className="columns is-desktop"
            style={{ margin: "15px", overflowX: "auto", whiteSpace: "nowrap" }}
          >
            {SMILES_LIST.map((smiles, index) => {
              // Dynamically calculate the width based on the molecule's size
              const dynamicWidth = subSMILES_LIST[index] ? subSMILES_LIST[index].length * 10 : moleculeWidth;
              return (
                <div className="column" key={smiles} style={{ width: `${dynamicWidth}px`, margin: "0 10px" }}>
                  <MoleculeStructure
                    id={smiles}
                    structure={smiles}
                    height={250}
                    width={moleculeWidth}
                    details={JSON.stringify({"atoms": subSMILES_LIST[index]})}
                    svgMode
                  />
                </div>
              );
            })}
          </div>
        </div>
      ) : null}
    </div>
  );
};

export default VisualOutput;

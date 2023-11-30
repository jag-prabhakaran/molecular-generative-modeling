import MoleculeStructure from "./MoleculeStructure";

const VisualOutput = ({ apiResponse }) => {
  const moleculeWidth = 300; // Set a default width for molecules
  const containerWidth = 90; // Set a default container width percentage

    console.log(apiResponse)
  return (
    <div>
      {apiResponse ? (
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

      {apiResponse ? (
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
            {apiResponse.map((responseObject, index) => {
              // Dynamically calculate the width based on the molecule's size
              const dynamicWidth = responseObject ? responseObject.length * 10 : moleculeWidth;
              var output_smile = null;
              var substructure = null;
              var logP = null;
              if(responseObject.hasOwnProperty('input_smile')){
                output_smile = responseObject.output_smile;
                substructure = responseObject.substructure_difference;  
                logP = responseObject.log_p;
              } else {
                output_smile = responseObject.smile;
                logP = responseObject.logP;
              }
              console.log(output_smile)
              console.log(substructure)
              console.log(logP)
              return (
                <div className="column" key={index} style={{ width: `${dynamicWidth}px`, margin: "0 10px" }}>
                  {substructure ? (
                  <MoleculeStructure
                    id={output_smile}
                    structure={output_smile}
                    height={250}
                    width={moleculeWidth}
                    details={JSON.stringify(
                      {
                        "atoms": substructure
                      })}
                    svgMode
                  />
                    ) : (
                  <MoleculeStructure
                    id={output_smile}
                    structure={output_smile}
                    height={250}
                    width={moleculeWidth}
                    details={"{}"}
                    svgMode
                  />
                  )
                    }
                </div>
              )
            }
            )}
          </div>
        </div>
      ) : null}
    </div>
  );
};

export default VisualOutput;

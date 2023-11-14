
const Properties = ({ logP, NumOfMolecules, onLogPChange, onNumOfMoleculesChange }) => {
    return (
        <div style={{
            position: "absolute",
            left: "75vw"
        }}>
            <div class="PropertiesText" style={{
                top: "-39vh",
                color: "black",
                fontSize: "2.5vh",
                fontFamily: "monospace",
                fontWeight: 400,
                wordWrap: "break-word",
                position: "relative"
            }}>
                Properties
            </div>
            <div class="logP" style={{
                background: "#8BCACC",
                borderRadius: "5px",
                width: "300px",
                height: "50px",
                top: "-37vh",
                position: "relative"
            }}>
                <div class="logPText" style={{
                    color: "black",
                    fontSize: "2.75vh",
                    fontFamily: "monospace",
                    fontWeight: 200,
                    wordWrap: "break-word",
                    top: "0.80vh",
                    left: "-2vw",
                    position: "relative"
                }}>
                    logP
                </div>
                <div class="logPValue">
                    <input style= {{
                        position: "relative",
                        top: "-3.5vh",
                        left: "4vw",
                        width: "50px",
                        height: "30px",
                        backgroundColor: "lightgray"
                    }}
                    type="number"
                    value={logP}
                    id="logPValue"
                    onChange={(e) => onLogPChange(e.target.value)}
                    placeholder="0"
                    >
                    </input>
                </div>
            </div>
            <div class="NumOfMolecules" style={{
                background: "#8BCACC",
                borderRadius: "5px",
                width: "300px",
                height: "50px",
                top: "-35vh",
                position: "relative"
            }}>
                <div class="NumOfMoleculesText" style={{
                    color: "black",
                    fontSize: "2vh",
                    fontFamily: "monospace",
                    fontWeight: 200,
                    wordWrap: "break-word",
                    top: "1.2vh",
                    left: "-2vw",
                    position: "relative"
                }}>
                    Num Of Molecules
                </div>
                <div class="NumOfMoleculesValue">
                    <input style= {{
                        position: "relative",
                        top: "-2.5vh",
                        left: "4vw",
                        width: "50px",
                        height: "30px",
                        backgroundColor: "lightgray"
                    }}
                    type="number"
                    value={NumOfMolecules}
                    id="NumOfMoleculesValue"
                    onChange={(e) => onNumOfMoleculesChange(e.target.value)}
                    placeholder="25"
                    >
                    </input>
                </div>
            </div>
        </div>
    )
}

export default Properties
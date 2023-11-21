
const Properties = ({ logPMax, onLogPMaxChange, logPMin, onLogPMinChange, NumOfMolecules, onNumOfMoleculesChange }) => {
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
                    fontSize: "2vh",
                    fontFamily: "monospace",
                    fontWeight: 200,
                    wordWrap: "break-word",
                    top: "1.5vh",
                    left: "1.5vw",
                    position: "absolute"
                }}>
                    logP
                </div>
                <div class="logPMaxValue">
                    <input style= {{
                        position: "absolute",
                        top: "0.7vh",
                        left: "5vw",
                        width: "50px",
                        height: "30px",
                        backgroundColor: "lightgray"
                    }}
                    type="number"
                    value={logPMax}
                    id="logPMax"
                    onChange={(e) => onLogPMaxChange(e.target.value)}>
                    </input>
                </div>
                <div class="ToText" style={{
                    color: "black",
                    fontSize: "2vh",
                    fontFamily: "monospace",
                    fontWeight: "150",
                    wordWrap: "break-word",
                    top: "1.5vh",
                    left: "9vw",
                    position: "absolute"
                }}>
                    to
                </div>
                <div class="logPMinValue">
                    <input style= {{
                        position: "absolute",
                        top: "0.7vh",
                        left: "11vw",
                        width: "50px",
                        height: "30px",
                        backgroundColor: "lightgray"
                    }}
                    type="number"
                    value={logPMin}
                    id="logPMin"
                    onChange={(e) => onLogPMinChange(e.target.value)}
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
                    >
                    </input>
                </div>
            </div>
        </div>
    )
}

export default Properties
import React from "react";

const ModelDescription = () => {
    return (
        <div>
            <div class = "model-description" style = {{
                position: 'absolute',
                width: '400px',
                height: '365px',
                top: '19vh',
                left: '2vw',
                backgroundColor: '#D0C1C1',
                borderRadius: '5px'
            }}>
                <div class = "model-des-text" style = {{
                    height: '50px',
                    width: '375px',
                    left: '0.25vw',
                    position: 'relative'
                }}>
                    <p style = {{
                        fontSize: '15px',
                        color: '#000',
                        fontFamily: "monospace"
                    }}>
                        This model employs a SMILES-based Recurrent Neural Network (RNN) generative model to 
                        achieve scaffold-constrained generation. 
                        The model takes as inputs :
                        <ul>
                            <li>a molecular scaffold in the form of a SMILES string with open positions marked by * (asterisks)
                                [there can be one or many open positions in a single scaffold SMILES string]</li>
                            <li>the upper bound of the number of molecules to be attempted to be generated. </li>
                        </ul>
                        The model outputs SMILES strings for the generated molecules.
                    </p>
                </div>
            </div>
        </div>
    );
};

export default ModelDescription
import React, { useEffect, useState } from 'react';


class MyComponent extends React.Component {
  constructor(props) {
    super(props);

    this.jsmeContainerRef = React.createRef();
    this.state = {
      // jsmeApplet: null, // Commented out for now
    };
  }

  componentDidMount() {
    // Initialize the JSME applet when the component mounts
    // const jsmeApplet = new JSApplet.JSME(this.jsmeContainerRef.current);
    // this.setState({ jsmeApplet }); // Commented out for now
  }

  showSMILES = () => {
    // const { jsmeApplet } = this.state;
    // if (jsmeApplet) {
    //   alert(jsmeApplet.smiles());
    // }
  }

  
  render() {
    return (
      <div>
        <head>
          <style>
            {`
              HTML {margin:0 !important; border:none !important;}
              .dragdrop-handle {cursor: move; user-select: none; -khtml-user-select: none; -moz-user-select: none;}
              .B {zoom: 1;}
              .dragdrop-dragging {zoom: normal;}
              .I {border: 1px dashed #1e90ff; margin: 0 !important; zoom: 1; z-index: 100;}
              .dragdrop-flow-panel-positioner {color: #1e90ff; display: inline; text-align: center; vertical-align: middle;}
              .dragdrop-proxy {background-color: #7af;}
              .dragdrop-selected, .dragdrop-dragging, .dragdrop-proxy {filter: alpha(opacity = 30); opacity: 0.3;}
              .dragdrop-movable-panel {z-index: 200; margin: 0 !important; border: none !important;}
            `}
          </style>
        </head>
        <body>
          <div className="desktop-main" style={{ width: '1440px', height: '1024px', position: 'relative', background: 'white' }}>
            <div className="ModelCheckpointext" style={{ color: 'black', fontSize: '30px', fontFamily: 'Inter', fontWeight: 400, top: '20px', left: '50px', wordWrap: 'break-word', position: 'absolute' }}>
              Model Checkpoint
            </div>
            <div className="Model Checkpoint" style={{ position: 'absolute', width: '400px', height: '365px', top: '160px', left: '50px', background: '#D0C1C1', borderRadius: '5px' }}>
              <div className="model-des-text" style={{ height: '50px', width: '375px', left: '15px', position: 'relative' }}>
                <p style={{ fontSize: '15px' }}>
                  This model employs a SMILES-based Recurrent Neural Network (RNN) generative model with a modified sampling procedure to achieve scaffold-constrained generation,
                  to help in the lead-optimization phase of drug discovery efforts. The model takes as inputs (1) a molecular scaffold in the form of a SMILES string with open
                  positions marked by * (asterisks) [there can be one or many open positions in a single scaffold SMILES string] and (2) the upper bound of the number of
                  molecules to be attempted to be generated. Note that the molecular generation process sometimes produces invalid SMILES strings and thus it often occurs
                  that the upper bound of the number of molecules to be generated is not reached (e.g.: when asking the model to generate 3000 molecules based on a scaffold,
                  around 2300-2500 candidate molecules are obtained). The model outputs SMILES strings for the generated molecules.
                </p>
              </div>
            </div>

            <div id="jsme_container" style={{ position: 'absolute', top: '5%', height: '20%', width: '45%', left: '38%' }}>
              <div tabIndex="0" style={{ padding: '0px', outline: '0px', width: '40%', height: '21.25%', position: 'absolute' }}></div>
            </div>
            {/*<div className="show-smiles" style={{ position: 'absolute', top: '430px', left: '550px' }}>
              <button type="button" onClick={() => alert(jsmeApplet.smiles())}>Show SMILES</button>
            </div>*/}
            <div className="Button-file-upload" style={{ position: 'absolute', top: '120px', height: '29px', width: '156px', left: '50px', background: '#8BCACC', borderRadius: '5px' }}>
              <button>Upload File</button>
            </div>
            <div className="generate button">
              <button id="generate" style={{ width: '243px', height: '52px', left: '550px', top: '475px', position: 'absolute', background: '#8BCACC', borderRadius: '5px' }}>
                GENERATE
              </button>
            </div>
            <div className="adme text" style={{ position: 'absolute', top: '50px', right: '-75px', color: 'black', fontSize: '30px', fontFamily: 'Inter', fontWeight: 400, wordWrap: 'break-word' }}>
              ADME Properties
            </div>
            <div className="logdinput" style={{ position: 'absolute', top: '100px', right: '-100px', width: '282px', height: '52px', background: '#8BCACC', borderRadius: '5px' }}></div>
            <div className="logdtext" style={{ position: 'absolute', top: '111.6px', right: '0px', width: '127.51px', height: '28.80px', color: 'black', fontSize: '24px', fontFamily: 'Inter', fontWeight: 400, wordWrap: 'break-word' }}>
              LOG D
            </div>
            <div className="logdvalue">
              <div>
                <input style={{ width: '54.77px', height: '29.60px', right: '-80px', top: '110.2px', position: 'absolute', background: '#D9D9D9', borderRadius: '2px' }} type="number" id="logd" value="%" />
                {/* Add your text or label here */}
                <span>Your Label or Text</span>
              </div>
            </div>
            <div className="qedinput" style={{ position: 'absolute', top: '175px', right: '-100px', width: '282px', height: '52px', background: '#8BCACC', borderRadius: '5px' }}></div>
            <div className="qedtext" style={{ position: 'absolute', top: '186px', right: '0px', width: '127.51px', height: '28.80px', color: 'black', fontSize: '24px', fontFamily: 'Inter', fontWeight: 400, wordWrap: 'break-word' }}>
              QED Score
            </div>
            <div className="visOutput-text" id="text" style={{ display: 'none', color: 'black', fontSize: '30px', fontFamily: 'Inter', fontWeight: 400, top: '550px', left: '550px', wordWrap: 'break-word', position: 'absolute' }}>
              Visual Output
            </div>
            <div className="visOutput" id="VO" style={{ display: 'none', position: 'absolute', top: '600px', left: '550px', background: '#D0C1C1', width: '680px', height: '300px', borderRadius: '5px' }}>
              <div className="visOutput-img" style={{ position: 'absolute', top: '20px', left: '50px' }}>
                <img className="output" src="molecule.png" style={{ height: '250px' }} />
              </div>
            </div>
          </div>
        </body>
      </div>
    );
  }
}

export default MyComponent;

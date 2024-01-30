"use client";
import React from "react";
import "miew/dist/Miew.min.css";
import { StandaloneStructServiceProvider } from "ketcher-standalone";
import { Editor } from "ketcher-react";
import "ketcher-react/dist/index.css";
// @ts-ignore
import Miew from "miew";

(window as any).Miew = Miew;

const structServiceProvider = new StandaloneStructServiceProvider();

export class KetcherComponent extends React.Component {
  ketcher: any;

  constructor(props: any) {
    super(props);
    this.handleOnInit = this.handleOnInit.bind(this);
  }

  handleOnInit = async (ketcher: any) => {
    this.ketcher = ketcher;
    (window as any).ketcher = ketcher;

    const initialData =
      "\n  Ketcher  9282116442D 1   1.00000     0.00000     0\n\n  6  6  0  0  0  0            999 V2000\n    9.5500  -11.7000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   10.4160  -12.2000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   10.4160  -13.2000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    9.5500  -13.7000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    8.6840  -13.2000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    8.6840  -12.2000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0     0  0\n  2  3  1  0     0  0\n  3  4  1  0     0  0\n  4  5  1  0     0  0\n  5  6  1  0     0  0\n  6  1  1  0     0  0\nM  END\n";
    this.ketcher.setMolecule(initialData);
  };

  render() {
    return (
      <Editor
        staticResourcesUrl={""}
        structServiceProvider={structServiceProvider}
        onInit={this.handleOnInit}
        errorHandler={(error: any) => console.error(error)}
         
      />
    );
  }
}

export default KetcherComponent;

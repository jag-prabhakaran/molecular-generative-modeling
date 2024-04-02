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
      "C1CCCCC1*";
    this.ketcher.setMolecule(initialData);
  };

  

  render() {
    if (typeof window === "undefined") {
      return <div>Server-side rendering is not supported</div>;
    } else {
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
}

export default KetcherComponent;

<!--- https://jojozhuang.github.io/tutorial/mermaid-cheat-sheet/ -->
<!--- https://mermaid-js.github.io/mermaid-live-editor/edit#pako:eNptkU1PwzAMhv9KlBMT6x-odkFskzjstNs0CbmJKVbzAfnQBF3_-9KyhNHhU_zYr_U67rmwEnnNhQLv1wStA300LMWTIQ2Krc5VxdZRdPd0S_79nh6wcfAH1-yRTGDQ4hzvgyPTshaNRHdbHCV-Bzo9HxazgoaAGU62J3v9D2BlaIPQPVtlXSn4E-ksTOlnBNHlfLidNy5W5lWjd0_f-GK2iKFgAWYD4V_99AW_hhprFSP_eiIlC3TRzLQ5-JJrdBpIprv0Q0rjh0w7byQF63j9BsrjkkMMdv9lBK-Di5ibrie8dg0XGq-Rag -->
<!--- Ctrl + K, V to split screen with output -->

```mermaid
flowchart LR

    __main__(__main__) --> A{Parse arguments}
    A --> ps
    A --> pj
    A --> pc

subgraph Split
    ps[process] --> split_file
    ps --> save_beam_lists
    ss[setup] --> gfs[get_files]
    ss --> gas[get_arc]
    split_file --> ss
    split_file --> split_ext
    split_file --> crop_file
    split_file --> update_beam_lists

end

subgraph Cross Correlation
    pc[process]

end

subgraph Join
    pj[process] --> join_file
    sj[setup] --> gfj[get_files]
    sj --> get_solutions
    sj --> gaj[get_arc]
    join_file --> sj
    join_file --> check_crop
    
end

```
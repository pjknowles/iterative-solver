<?xml version='1.0' encoding='UTF-8' standalone='yes' ?>
<tagfile>
  <compound kind="struct">
    <name>molpro::profiler::detail::AccessCalls</name>
    <filename>structmolpro_1_1profiler_1_1detail_1_1AccessCalls.html</filename>
    <member kind="function">
      <type>double</type>
      <name>operator()</name>
      <anchorfile>structmolpro_1_1profiler_1_1detail_1_1AccessCalls.html</anchorfile>
      <anchor>a0e1d9458d1d8e25ef508498ca1e1d04c</anchor>
      <arglist>(const TreePath &amp;t)</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>molpro::profiler::detail::AccessCPU</name>
    <filename>structmolpro_1_1profiler_1_1detail_1_1AccessCPU.html</filename>
    <member kind="function">
      <type>double</type>
      <name>operator()</name>
      <anchorfile>structmolpro_1_1profiler_1_1detail_1_1AccessCPU.html</anchorfile>
      <anchor>a5bd501ceac97b99963f876daefdc8226</anchor>
      <arglist>(const TreePath &amp;t)</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>molpro::profiler::detail::AccessOperations</name>
    <filename>structmolpro_1_1profiler_1_1detail_1_1AccessOperations.html</filename>
    <member kind="function">
      <type>double</type>
      <name>operator()</name>
      <anchorfile>structmolpro_1_1profiler_1_1detail_1_1AccessOperations.html</anchorfile>
      <anchor>af135568aa0e88453beefea0de64f077e</anchor>
      <arglist>(const TreePath &amp;t)</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>molpro::profiler::detail::AccessWall</name>
    <filename>structmolpro_1_1profiler_1_1detail_1_1AccessWall.html</filename>
    <member kind="function">
      <type>double</type>
      <name>operator()</name>
      <anchorfile>structmolpro_1_1profiler_1_1detail_1_1AccessWall.html</anchorfile>
      <anchor>a780881aebc4449dee7f3329bf59e6269</anchor>
      <arglist>(const TreePath &amp;t)</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>molpro::profiler::detail::Compare</name>
    <filename>structmolpro_1_1profiler_1_1detail_1_1Compare.html</filename>
    <templarg></templarg>
    <member kind="function">
      <type>bool</type>
      <name>operator()</name>
      <anchorfile>structmolpro_1_1profiler_1_1detail_1_1Compare.html</anchorfile>
      <anchor>a7b20b0e27a091b98f54acfc5b5da4e50</anchor>
      <arglist>(const TreePath &amp;l, const TreePath &amp;r) const</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>molpro::profiler::Counter</name>
    <filename>classmolpro_1_1profiler_1_1Counter.html</filename>
    <member kind="function">
      <type></type>
      <name>Counter</name>
      <anchorfile>classmolpro_1_1profiler_1_1Counter.html</anchorfile>
      <anchor>aa04f033daefaf7d484734f176ba79e1b</anchor>
      <arglist>()=default</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Counter</name>
      <anchorfile>classmolpro_1_1profiler_1_1Counter.html</anchorfile>
      <anchor>ac6dd698ff96df251d4191f4d1094d059</anchor>
      <arglist>(bool with_cpu_time, bool with_wall_time)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Counter</name>
      <anchorfile>classmolpro_1_1profiler_1_1Counter.html</anchorfile>
      <anchor>a3fec52ce9a1583b3da21098e46d37df2</anchor>
      <arglist>(size_t call_count_, size_t operation_count_, double wall_time_, double cpu_time_, bool with_cpu_time, bool with_wall_time)</arglist>
    </member>
    <member kind="function">
      <type>Counter &amp;</type>
      <name>start</name>
      <anchorfile>classmolpro_1_1profiler_1_1Counter.html</anchorfile>
      <anchor>a045bab06889c88d495df04e2f9549d06</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Counter &amp;</type>
      <name>stop</name>
      <anchorfile>classmolpro_1_1profiler_1_1Counter.html</anchorfile>
      <anchor>ac282bf4bfd1bd0f9b8e3e7bf49aec78b</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Counter &amp;</type>
      <name>reset</name>
      <anchorfile>classmolpro_1_1profiler_1_1Counter.html</anchorfile>
      <anchor>a0d0e79d1329c3bf5c9b0d1de71ec7cec</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>add_operations</name>
      <anchorfile>classmolpro_1_1profiler_1_1Counter.html</anchorfile>
      <anchor>a1c45213b0bc43c474c736180debb095d</anchor>
      <arglist>(size_t ops)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator+=</name>
      <anchorfile>classmolpro_1_1profiler_1_1Counter.html</anchorfile>
      <anchor>a2ec28e4a7efe5ddec22ade05d8908490</anchor>
      <arglist>(const Counter &amp;other)</arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>get_call_count</name>
      <anchorfile>classmolpro_1_1profiler_1_1Counter.html</anchorfile>
      <anchor>a7dd38b22b6c79eef610b672f56aeb102</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>get_operation_count</name>
      <anchorfile>classmolpro_1_1profiler_1_1Counter.html</anchorfile>
      <anchor>a167237d369446b5ae4b485b30cfab6b6</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>const Timer &amp;</type>
      <name>get_cpu</name>
      <anchorfile>classmolpro_1_1profiler_1_1Counter.html</anchorfile>
      <anchor>a7cdc226239439830264bcaaff07fb7a8</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>const Timer &amp;</type>
      <name>get_wall</name>
      <anchorfile>classmolpro_1_1profiler_1_1Counter.html</anchorfile>
      <anchor>adbf2582f654d7bdcb1dcd7cae25ea8c1</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>size_t</type>
      <name>call_count</name>
      <anchorfile>classmolpro_1_1profiler_1_1Counter.html</anchorfile>
      <anchor>aa9f4ab9e1c9e410f673ddbc081fd4783</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>size_t</type>
      <name>operation_count</name>
      <anchorfile>classmolpro_1_1profiler_1_1Counter.html</anchorfile>
      <anchor>ae417f08433018267b65ddebe77852719</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Timer</type>
      <name>cpu</name>
      <anchorfile>classmolpro_1_1profiler_1_1Counter.html</anchorfile>
      <anchor>afbaa26e1b0f2f6ba8009a9ed5923f245</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Timer</type>
      <name>wall</name>
      <anchorfile>classmolpro_1_1profiler_1_1Counter.html</anchorfile>
      <anchor>a99620a5bacd3e931858d4656c50aecf7</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>molpro::profiler::dotgraph::GraphEntry</name>
    <filename>classmolpro_1_1profiler_1_1dotgraph_1_1GraphEntry.html</filename>
    <member kind="function">
      <type></type>
      <name>GraphEntry</name>
      <anchorfile>classmolpro_1_1profiler_1_1dotgraph_1_1GraphEntry.html</anchorfile>
      <anchor>a9b52746443308c8b8a7ff81212339faf</anchor>
      <arglist>(EntryType entry_type, std::string name, double runtime, int calls, double total_time, int operations=-1, std::string name_to=&quot;&quot;)</arglist>
    </member>
    <member kind="function">
      <type>std::pair&lt; std::string, std::string &gt;</type>
      <name>get_colours</name>
      <anchorfile>classmolpro_1_1profiler_1_1dotgraph_1_1GraphEntry.html</anchorfile>
      <anchor>a27009f70bf7c7740e01fe9c992a92f01</anchor>
      <arglist>(int hot[3], int cool[3], double total_time)</arglist>
    </member>
    <member kind="variable">
      <type>EntryType</type>
      <name>entry_type</name>
      <anchorfile>classmolpro_1_1profiler_1_1dotgraph_1_1GraphEntry.html</anchorfile>
      <anchor>acc5e3f71df0aa85f6f8268fc3a3acca5</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::string</type>
      <name>name</name>
      <anchorfile>classmolpro_1_1profiler_1_1dotgraph_1_1GraphEntry.html</anchorfile>
      <anchor>a7f0aef10ddc7d331b155aade1f0e84d7</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>runtime</name>
      <anchorfile>classmolpro_1_1profiler_1_1dotgraph_1_1GraphEntry.html</anchorfile>
      <anchor>a6d8caabf8b7aca287bcbb40169d25b27</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>calls</name>
      <anchorfile>classmolpro_1_1profiler_1_1dotgraph_1_1GraphEntry.html</anchorfile>
      <anchor>a45b8d66eb8b49c4aea7a0ef35493b963</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::string</type>
      <name>name_to</name>
      <anchorfile>classmolpro_1_1profiler_1_1dotgraph_1_1GraphEntry.html</anchorfile>
      <anchor>a9248fd2eb409b7636ec30727f3e47205</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::string</type>
      <name>fontcolour</name>
      <anchorfile>classmolpro_1_1profiler_1_1dotgraph_1_1GraphEntry.html</anchorfile>
      <anchor>a0fdd82c75a18fd746d6d5873284d4aa2</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>operations</name>
      <anchorfile>classmolpro_1_1profiler_1_1dotgraph_1_1GraphEntry.html</anchorfile>
      <anchor>a0f7c348d94630a806490824ca320c21a</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>molpro::profiler::Node</name>
    <filename>classmolpro_1_1profiler_1_1Node.html</filename>
    <templarg></templarg>
    <member kind="function">
      <type></type>
      <name>~Node</name>
      <anchorfile>classmolpro_1_1profiler_1_1Node.html</anchorfile>
      <anchor>aff949cced22954074f78bfe430d84250</anchor>
      <arglist>()=default</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Node</name>
      <anchorfile>classmolpro_1_1profiler_1_1Node.html</anchorfile>
      <anchor>a9fd8e958e669ccaf3fff71a4e24c8baa</anchor>
      <arglist>(Node&lt; Counter &gt; &amp;&amp;) noexcept=default</arglist>
    </member>
    <member kind="function">
      <type>Node&lt; Counter &gt; &amp;</type>
      <name>operator=</name>
      <anchorfile>classmolpro_1_1profiler_1_1Node.html</anchorfile>
      <anchor>a344bac5389fb9f12438d3df02a24edce</anchor>
      <arglist>(Node&lt; Counter &gt; &amp;&amp;) noexcept=default</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Node</name>
      <anchorfile>classmolpro_1_1profiler_1_1Node.html</anchorfile>
      <anchor>a1971a26b3c38e44c8228f7e1f7d8daf9</anchor>
      <arglist>(const Node&lt; Counter &gt; &amp;)=delete</arglist>
    </member>
    <member kind="function">
      <type>Node&lt; Counter &gt; &amp;</type>
      <name>operator=</name>
      <anchorfile>classmolpro_1_1profiler_1_1Node.html</anchorfile>
      <anchor>a20731d0b1510610a88cd22a999b029dd</anchor>
      <arglist>(const Node&lt; Counter &gt; &amp;)=delete</arglist>
    </member>
    <member kind="function">
      <type>std::shared_ptr&lt; Node&lt; Counter &gt; &gt;</type>
      <name>walk</name>
      <anchorfile>classmolpro_1_1profiler_1_1Node.html</anchorfile>
      <anchor>a19f49c5027c223af3e229c50c592b261</anchor>
      <arglist>(ForwardIt start_name, ForwardIt end_name)</arglist>
    </member>
    <member kind="function">
      <type>std::shared_ptr&lt; Node&lt; Counter &gt; &gt;</type>
      <name>walk</name>
      <anchorfile>classmolpro_1_1profiler_1_1Node.html</anchorfile>
      <anchor>a0c14484f400e493fc13b14d96c06c12f</anchor>
      <arglist>(const std::list&lt; std::string &gt; &amp;path_to_node)</arglist>
    </member>
    <member kind="function">
      <type>std::shared_ptr&lt; Node&lt; Counter &gt; &gt;</type>
      <name>child</name>
      <anchorfile>classmolpro_1_1profiler_1_1Node.html</anchorfile>
      <anchor>a065f0e3d02f1e9be69df873391edae54</anchor>
      <arglist>(const std::string &amp;child_name)</arglist>
    </member>
    <member kind="function">
      <type>std::shared_ptr&lt; Node&lt; Counter &gt; &gt;</type>
      <name>find_parent</name>
      <anchorfile>classmolpro_1_1profiler_1_1Node.html</anchorfile>
      <anchor>a7d40f6320d59b8a3c857a855cd83c19f</anchor>
      <arglist>(const std::string &amp;parent_name)</arglist>
    </member>
    <member kind="function">
      <type>std::shared_ptr&lt; Node&lt; Counter &gt; &gt;</type>
      <name>walk_up</name>
      <anchorfile>classmolpro_1_1profiler_1_1Node.html</anchorfile>
      <anchor>ad21e99b5ef2e90fa9fd461c0489a0114</anchor>
      <arglist>(int n)</arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>count_nodes</name>
      <anchorfile>classmolpro_1_1profiler_1_1Node.html</anchorfile>
      <anchor>acd0c3a6cd57c281dbc666a3db1401f4a</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static std::shared_ptr&lt; Node&lt; Counter &gt; &gt;</type>
      <name>make_root</name>
      <anchorfile>classmolpro_1_1profiler_1_1Node.html</anchorfile>
      <anchor>a607c98e74e445367501cbbba0b614d01</anchor>
      <arglist>(const std::string &amp;name, const Counter &amp;counter)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static std::shared_ptr&lt; Node&lt; Counter &gt; &gt;</type>
      <name>add_child</name>
      <anchorfile>classmolpro_1_1profiler_1_1Node.html</anchorfile>
      <anchor>ae15812662da535d34925cb5d84bfb6e6</anchor>
      <arglist>(const std::string &amp;child_name, const Counter &amp;child_counter, const std::shared_ptr&lt; Node&lt; Counter &gt;&gt; &amp;parent)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static std::shared_ptr&lt; Node&lt; Counter &gt; &gt;</type>
      <name>deep_copy</name>
      <anchorfile>classmolpro_1_1profiler_1_1Node.html</anchorfile>
      <anchor>a172c8ef6d50161afd3d605aff3e2d49b</anchor>
      <arglist>(const std::shared_ptr&lt; Node&lt; Counter &gt;&gt; &amp;subtree, std::shared_ptr&lt; Node&lt; Counter &gt;&gt; parent)</arglist>
    </member>
    <member kind="variable">
      <type>std::string</type>
      <name>name</name>
      <anchorfile>classmolpro_1_1profiler_1_1Node.html</anchorfile>
      <anchor>a6df8b209b0dc6b0b71c21dafaf2d83bf</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>Counter</type>
      <name>counter</name>
      <anchorfile>classmolpro_1_1profiler_1_1Node.html</anchorfile>
      <anchor>a4409a690d20fa5224d1e4ec656be8f4e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::shared_ptr&lt; Node&lt; Counter &gt; &gt;</type>
      <name>parent</name>
      <anchorfile>classmolpro_1_1profiler_1_1Node.html</anchorfile>
      <anchor>a87140036a8602a54fab125a88f7b32b0</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::map&lt; std::string, std::shared_ptr&lt; Node&lt; Counter &gt; &gt; &gt;</type>
      <name>children</name>
      <anchorfile>classmolpro_1_1profiler_1_1Node.html</anchorfile>
      <anchor>a5dbdef08dfb1a76a51e0ab35890926b8</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>molpro::profiler::NodePathError</name>
    <filename>structmolpro_1_1profiler_1_1NodePathError.html</filename>
  </compound>
  <compound kind="struct">
    <name>molpro::profiler::detail::None</name>
    <filename>structmolpro_1_1profiler_1_1detail_1_1None.html</filename>
    <member kind="function">
      <type>double</type>
      <name>operator()</name>
      <anchorfile>structmolpro_1_1profiler_1_1detail_1_1None.html</anchorfile>
      <anchor>aa45256fef921ad5f3f7712ecd78a4c6e</anchor>
      <arglist>(const TreePath &amp;t)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>operator()</name>
      <anchorfile>structmolpro_1_1profiler_1_1detail_1_1None.html</anchorfile>
      <anchor>aa1a693d3482687c05a1fb7b05a0de804</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="interface">
    <name>profilerf::profiler</name>
    <filename>structprofilerf_1_1profiler.html</filename>
    <member kind="function">
      <type>PROCEDURE</type>
      <name>start</name>
      <anchorfile>structprofilerf_1_1profiler.html</anchorfile>
      <anchor>a3e0782a7c0064f2ab915da6c0c743235</anchor>
      <arglist>=&gt; profilerstartf</arglist>
    </member>
    <member kind="function">
      <type>PROCEDURE</type>
      <name>stop</name>
      <anchorfile>structprofilerf_1_1profiler.html</anchorfile>
      <anchor>a58e3ce95e466c417bd148d32f35017d4</anchor>
      <arglist>=&gt; profilerstopf</arglist>
    </member>
    <member kind="function">
      <type>PROCEDURE</type>
      <name>active</name>
      <anchorfile>structprofilerf_1_1profiler.html</anchorfile>
      <anchor>aa8c588e2033df32b86a5141bf6d6e437</anchor>
      <arglist>=&gt; profileractivef</arglist>
    </member>
    <member kind="function">
      <type>PROCEDURE</type>
      <name>print</name>
      <anchorfile>structprofilerf_1_1profiler.html</anchorfile>
      <anchor>aadde1e394375eae47707367fb4ca1296</anchor>
      <arglist>=&gt; profilerprintf</arglist>
    </member>
    <member kind="function">
      <type>PROCEDURE</type>
      <name>destroy</name>
      <anchorfile>structprofilerf_1_1profiler.html</anchorfile>
      <anchor>a776c8a7adb921fb73699ec5141af37dd</anchor>
      <arglist>=&gt; profilerdestroyf</arglist>
    </member>
    <member kind="function">
      <type>PROCEDURE</type>
      <name>dotgraph</name>
      <anchorfile>structprofilerf_1_1profiler.html</anchorfile>
      <anchor>ae1d67163b39c8d0a2bb9eb682d708dc1</anchor>
      <arglist>=&gt; profilerdotgraphf</arglist>
    </member>
    <member kind="function">
      <type>FINAL</type>
      <name>destructor</name>
      <anchorfile>structprofilerf_1_1profiler.html</anchorfile>
      <anchor>a75d48c3b53ca9ea23b63ee2666c562ca</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>type(profiler) function</type>
      <name>profilernewf</name>
      <anchorfile>structprofilerf_1_1profiler.html</anchorfile>
      <anchor>a884ba4d385ef346cc3da0e44f7159e5e</anchor>
      <arglist>(name, sort, level, comm, cpu)</arglist>
    </member>
    <member kind="variable">
      <type>type(c_ptr)</type>
      <name>handle</name>
      <anchorfile>structprofilerf_1_1profiler.html</anchorfile>
      <anchor>ab0fc2e7b48fa188b792f578e863e2afd</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>molpro::profiler::Profiler</name>
    <filename>classmolpro_1_1profiler_1_1Profiler.html</filename>
    <class kind="struct">molpro::profiler::Profiler::Proxy</class>
    <member kind="function">
      <type></type>
      <name>Profiler</name>
      <anchorfile>classmolpro_1_1profiler_1_1Profiler.html</anchorfile>
      <anchor>a506158e05edf11a9d734a0fe1097a422</anchor>
      <arglist>(std::string description_, bool with_wall=true, bool with_cpu=false)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Profiler</name>
      <anchorfile>classmolpro_1_1profiler_1_1Profiler.html</anchorfile>
      <anchor>aec32d1f8dd5a31526b3c397c7c9304d5</anchor>
      <arglist>(Profiler &amp;&amp;)=default</arglist>
    </member>
    <member kind="function">
      <type>Profiler &amp;</type>
      <name>operator=</name>
      <anchorfile>classmolpro_1_1profiler_1_1Profiler.html</anchorfile>
      <anchor>a5aa955177d5c66b59ef5c7a655d39cec</anchor>
      <arglist>(Profiler &amp;&amp;)=default</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~Profiler</name>
      <anchorfile>classmolpro_1_1profiler_1_1Profiler.html</anchorfile>
      <anchor>a9a833718c09eaf55066af3d8a2c23f13</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Profiler</name>
      <anchorfile>classmolpro_1_1profiler_1_1Profiler.html</anchorfile>
      <anchor>a967e7a44b23950697443038719539e9b</anchor>
      <arglist>()=delete</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Profiler</name>
      <anchorfile>classmolpro_1_1profiler_1_1Profiler.html</anchorfile>
      <anchor>ae7404abd6d208a1c3cd0f1ffe466cf9c</anchor>
      <arglist>(const Profiler &amp;)=delete</arglist>
    </member>
    <member kind="function">
      <type>Profiler &amp;</type>
      <name>operator=</name>
      <anchorfile>classmolpro_1_1profiler_1_1Profiler.html</anchorfile>
      <anchor>a92b04f4d8fb586675c205025c7af71b9</anchor>
      <arglist>(const Profiler &amp;)=delete</arglist>
    </member>
    <member kind="function">
      <type>const std::string &amp;</type>
      <name>description</name>
      <anchorfile>classmolpro_1_1profiler_1_1Profiler.html</anchorfile>
      <anchor>a891098d9b1402b5ae72554a1ce68bd9e</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>get_max_depth</name>
      <anchorfile>classmolpro_1_1profiler_1_1Profiler.html</anchorfile>
      <anchor>a6d9491d3523922d3b262f3b3d807f96e</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_max_depth</name>
      <anchorfile>classmolpro_1_1profiler_1_1Profiler.html</anchorfile>
      <anchor>a555c0251673b04cbc2e5ee814d23cda9</anchor>
      <arglist>(int new_max_depth)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>get_current_depth</name>
      <anchorfile>classmolpro_1_1profiler_1_1Profiler.html</anchorfile>
      <anchor>ae06994fb3ee65581e41fb9b3424a196f</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>Profiler &amp;</type>
      <name>start</name>
      <anchorfile>classmolpro_1_1profiler_1_1Profiler.html</anchorfile>
      <anchor>a73fc3106d01b5766bb250b05cc2d39e8</anchor>
      <arglist>(const std::string &amp;name)</arglist>
    </member>
    <member kind="function">
      <type>Profiler &amp;</type>
      <name>stop</name>
      <anchorfile>classmolpro_1_1profiler_1_1Profiler.html</anchorfile>
      <anchor>abc8010ddc06df088ac036b0ed2697737</anchor>
      <arglist>(const std::string &amp;name=&quot;&quot;)</arglist>
    </member>
    <member kind="function">
      <type>Profiler &amp;</type>
      <name>stop_all</name>
      <anchorfile>classmolpro_1_1profiler_1_1Profiler.html</anchorfile>
      <anchor>a5160bc6c72dc1f9184944ee837755618</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>profiler::Counter &amp;</type>
      <name>counter</name>
      <anchorfile>classmolpro_1_1profiler_1_1Profiler.html</anchorfile>
      <anchor>a2c2630dc28710f397982f7f55c3e67c0</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Profiler &amp;</type>
      <name>reset</name>
      <anchorfile>classmolpro_1_1profiler_1_1Profiler.html</anchorfile>
      <anchor>ac81d8f1f81d48e46105b5d6b2a46b1f3</anchor>
      <arglist>(const std::string &amp;name)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator+=</name>
      <anchorfile>classmolpro_1_1profiler_1_1Profiler.html</anchorfile>
      <anchor>a65992a880172726fde788d2cfd468fde</anchor>
      <arglist>(size_t operations)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator++</name>
      <anchorfile>classmolpro_1_1profiler_1_1Profiler.html</anchorfile>
      <anchor>a2864d7dfd267777b37f5437706f0b579</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>operator++</name>
      <anchorfile>classmolpro_1_1profiler_1_1Profiler.html</anchorfile>
      <anchor>a4dcc9cadba872475bbbbf3792eacca5f</anchor>
      <arglist>(int)</arglist>
    </member>
    <member kind="function">
      <type>Proxy</type>
      <name>push</name>
      <anchorfile>classmolpro_1_1profiler_1_1Profiler.html</anchorfile>
      <anchor>ac38692f6d56f4ea12acebe8195b1b20d</anchor>
      <arglist>(const std::string &amp;name)</arglist>
    </member>
    <member kind="function">
      <type>std::string</type>
      <name>str</name>
      <anchorfile>classmolpro_1_1profiler_1_1Profiler.html</anchorfile>
      <anchor>a87abba331dc948e200de310fd40d78e7</anchor>
      <arglist>(bool cumulative=true, profiler::SortBy sort_by=profiler::SortBy::wall) const</arglist>
    </member>
    <member kind="function">
      <type>std::string</type>
      <name>dotgraph</name>
      <anchorfile>classmolpro_1_1profiler_1_1Profiler.html</anchorfile>
      <anchor>a2b896c4cbc5312db8ea4d5f88008063e</anchor>
      <arglist>(std::string path, double threshold=0.01, bool cumulative=true, int hot[3]=hot_default, int cool[3]=cool_default, SortBy sort_by=profiler::SortBy::none, std::vector&lt; std::pair&lt; double, double &gt;&gt; heat_adjust=default_heat_adjust, bool get_percentage_time=false)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static std::shared_ptr&lt; Profiler &gt;</type>
      <name>single</name>
      <anchorfile>classmolpro_1_1profiler_1_1Profiler.html</anchorfile>
      <anchor>ae39b1e0c136785340f79c4a91897b001</anchor>
      <arglist>(const std::string &amp;description_, bool with_wall=true, bool with_cpu=false)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static std::shared_ptr&lt; Profiler &gt;</type>
      <name>single</name>
      <anchorfile>classmolpro_1_1profiler_1_1Profiler.html</anchorfile>
      <anchor>aea6108170dddaa5caafb2df4ec42db5d</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>erase</name>
      <anchorfile>classmolpro_1_1profiler_1_1Profiler.html</anchorfile>
      <anchor>a337a11c110c072b8c2a778dbad651be7</anchor>
      <arglist>(const std::string &amp;description)</arglist>
    </member>
    <member kind="variable">
      <type>std::shared_ptr&lt; profiler::Node&lt; profiler::Counter &gt; &gt;</type>
      <name>root</name>
      <anchorfile>classmolpro_1_1profiler_1_1Profiler.html</anchorfile>
      <anchor>a0ab99cd9b2d8b922b6529d19a42c30c5</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::shared_ptr&lt; profiler::Node&lt; profiler::Counter &gt; &gt;</type>
      <name>active_node</name>
      <anchorfile>classmolpro_1_1profiler_1_1Profiler.html</anchorfile>
      <anchor>a4e780f947690671a2948030b697f4e0d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>std::string</type>
      <name>m_description</name>
      <anchorfile>classmolpro_1_1profiler_1_1Profiler.html</anchorfile>
      <anchor>a4e8d9cdbc81cf2c060933670f1726e81</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>std::string</type>
      <name>m_root_name</name>
      <anchorfile>classmolpro_1_1profiler_1_1Profiler.html</anchorfile>
      <anchor>a545641f4683b65386f8688ede8e61392</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>int</type>
      <name>m_max_depth</name>
      <anchorfile>classmolpro_1_1profiler_1_1Profiler.html</anchorfile>
      <anchor>a39d52ee3307a2940e636ffcb5c3cf34e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>int</type>
      <name>m_current_depth</name>
      <anchorfile>classmolpro_1_1profiler_1_1Profiler.html</anchorfile>
      <anchor>abda688ae7d8d0050bf1588bbebfaa84d</anchor>
      <arglist></arglist>
    </member>
    <member kind="friend">
      <type>friend std::ostream &amp;</type>
      <name>operator&lt;&lt;</name>
      <anchorfile>classmolpro_1_1profiler_1_1Profiler.html</anchorfile>
      <anchor>aaffc7620e7e7ed86567cfbf30d686cea</anchor>
      <arglist>(std::ostream &amp;os, const Profiler &amp;obj)</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>molpro::profiler::Profiler::Proxy</name>
    <filename>structmolpro_1_1profiler_1_1Profiler_1_1Proxy.html</filename>
    <member kind="function">
      <type></type>
      <name>Proxy</name>
      <anchorfile>structmolpro_1_1profiler_1_1Profiler_1_1Proxy.html</anchorfile>
      <anchor>ac6fedc2d3892135bc89c1688015652aa</anchor>
      <arglist>(Profiler &amp;prof, const std::string &amp;name)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~Proxy</name>
      <anchorfile>structmolpro_1_1profiler_1_1Profiler_1_1Proxy.html</anchorfile>
      <anchor>aabc2bebffdde9361c0424c9f2860620c</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Proxy</type>
      <name>push</name>
      <anchorfile>structmolpro_1_1profiler_1_1Profiler_1_1Proxy.html</anchorfile>
      <anchor>af736944486dafc9a533f823ba656e957</anchor>
      <arglist>(const std::string &amp;name)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator+=</name>
      <anchorfile>structmolpro_1_1profiler_1_1Profiler_1_1Proxy.html</anchorfile>
      <anchor>ab1667157c570954bf1ffe74c7ffc2aa7</anchor>
      <arglist>(size_t operations)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator++</name>
      <anchorfile>structmolpro_1_1profiler_1_1Profiler_1_1Proxy.html</anchorfile>
      <anchor>ae135f9250485e847d060b51bf628790f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>operator++</name>
      <anchorfile>structmolpro_1_1profiler_1_1Profiler_1_1Proxy.html</anchorfile>
      <anchor>a9a3e6ff08c1fcae3738c0c6a45c41256</anchor>
      <arglist>(int)</arglist>
    </member>
    <member kind="variable">
      <type>Profiler &amp;</type>
      <name>prof</name>
      <anchorfile>structmolpro_1_1profiler_1_1Profiler_1_1Proxy.html</anchorfile>
      <anchor>a070687ea3da0dbc7d8f419132aa57d78</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>molpro::profiler::Timer</name>
    <filename>classmolpro_1_1profiler_1_1Timer.html</filename>
    <member kind="enumeration">
      <type></type>
      <name>Type</name>
      <anchorfile>classmolpro_1_1profiler_1_1Timer.html</anchorfile>
      <anchor>afba52e6cbab6bd528476f2f43d98fda8</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>cpu</name>
      <anchorfile>classmolpro_1_1profiler_1_1Timer.html</anchorfile>
      <anchor>afba52e6cbab6bd528476f2f43d98fda8a41e25934fb914bf1c819356146afd160</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>wall</name>
      <anchorfile>classmolpro_1_1profiler_1_1Timer.html</anchorfile>
      <anchor>afba52e6cbab6bd528476f2f43d98fda8aba4095fd99d7a4aa7c3ce9a377a7d6ce</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>cpu</name>
      <anchorfile>classmolpro_1_1profiler_1_1Timer.html</anchorfile>
      <anchor>afba52e6cbab6bd528476f2f43d98fda8a41e25934fb914bf1c819356146afd160</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>wall</name>
      <anchorfile>classmolpro_1_1profiler_1_1Timer.html</anchorfile>
      <anchor>afba52e6cbab6bd528476f2f43d98fda8aba4095fd99d7a4aa7c3ce9a377a7d6ce</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Timer</name>
      <anchorfile>classmolpro_1_1profiler_1_1Timer.html</anchorfile>
      <anchor>a9d9f382f87bd64d9189e5ce3965399ea</anchor>
      <arglist>(Type type, bool is_dummy)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Timer</name>
      <anchorfile>classmolpro_1_1profiler_1_1Timer.html</anchorfile>
      <anchor>a2b0c4272fe373b7abaee902b338adb02</anchor>
      <arglist>(double cumulative_time, Type type, bool is_dummy)</arglist>
    </member>
    <member kind="function">
      <type>Timer &amp;</type>
      <name>start</name>
      <anchorfile>classmolpro_1_1profiler_1_1Timer.html</anchorfile>
      <anchor>a6baa7c7178af4c55fe6adeaace7740b1</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Timer &amp;</type>
      <name>stop</name>
      <anchorfile>classmolpro_1_1profiler_1_1Timer.html</anchorfile>
      <anchor>ac72431aa164dec83bc81a8938da13ffd</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator+=</name>
      <anchorfile>classmolpro_1_1profiler_1_1Timer.html</anchorfile>
      <anchor>a85ca379be4a1ddd13fe694c1c8c950b2</anchor>
      <arglist>(const Timer &amp;other)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>start_time</name>
      <anchorfile>classmolpro_1_1profiler_1_1Timer.html</anchorfile>
      <anchor>ae39de3b76da661d8d3f6a715e7a7adca</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>stop_time</name>
      <anchorfile>classmolpro_1_1profiler_1_1Timer.html</anchorfile>
      <anchor>a236204959708a63d1cd501218e0a6671</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>cumulative_time</name>
      <anchorfile>classmolpro_1_1profiler_1_1Timer.html</anchorfile>
      <anchor>a361230d93a3aafd78998828e81653fad</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>stopped</name>
      <anchorfile>classmolpro_1_1profiler_1_1Timer.html</anchorfile>
      <anchor>a9de188061be00a37ac08b308f3aa4674</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>dummy</name>
      <anchorfile>classmolpro_1_1profiler_1_1Timer.html</anchorfile>
      <anchor>ad0f34fd523173a272c8ed7e20f8083f7</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>reset</name>
      <anchorfile>classmolpro_1_1profiler_1_1Timer.html</anchorfile>
      <anchor>a9bd2b1093d392ae88740485d094ea45c</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Type</type>
      <name>type</name>
      <anchorfile>classmolpro_1_1profiler_1_1Timer.html</anchorfile>
      <anchor>ab61998075350f0b424f3c0abcb49437c</anchor>
      <arglist>() const</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>molpro::profiler::detail::TreePath</name>
    <filename>structmolpro_1_1profiler_1_1detail_1_1TreePath.html</filename>
    <member kind="function">
      <type></type>
      <name>TreePath</name>
      <anchorfile>structmolpro_1_1profiler_1_1detail_1_1TreePath.html</anchorfile>
      <anchor>ada8d4bacd52b41d790896b4bd460b255</anchor>
      <arglist>(std::shared_ptr&lt; Node&lt; Counter &gt;&gt; node, bool cumulative)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static std::list&lt; TreePath &gt;</type>
      <name>convert_tree_to_paths</name>
      <anchorfile>structmolpro_1_1profiler_1_1detail_1_1TreePath.html</anchorfile>
      <anchor>a9009dec346015ab8ac7e40ff2366338d</anchor>
      <arglist>(const std::shared_ptr&lt; Node&lt; Counter &gt;&gt; &amp;root, bool cumulative, SortBy sort_by)</arglist>
    </member>
    <member kind="variable">
      <type>Counter</type>
      <name>counter</name>
      <anchorfile>structmolpro_1_1profiler_1_1detail_1_1TreePath.html</anchorfile>
      <anchor>aababc1fcd2ff5488f226399f4c1ef050</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::list&lt; std::string &gt;</type>
      <name>path</name>
      <anchorfile>structmolpro_1_1profiler_1_1detail_1_1TreePath.html</anchorfile>
      <anchor>aed5c02639d52ec44aa114b973fa04e98</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>size_t</type>
      <name>depth</name>
      <anchorfile>structmolpro_1_1profiler_1_1detail_1_1TreePath.html</anchorfile>
      <anchor>a2b542c4f02b397509df19d6e8d6a3a1b</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>molpro::profiler::WeakSingleton</name>
    <filename>structmolpro_1_1profiler_1_1WeakSingleton.html</filename>
    <templarg></templarg>
    <member kind="typedef">
      <type>std::tuple&lt; std::string, std::weak_ptr&lt; Object &gt;, Object * &gt;</type>
      <name>key_t</name>
      <anchorfile>structmolpro_1_1profiler_1_1WeakSingleton.html</anchorfile>
      <anchor>a914ac8e10632f06314f5bf81e4845a58</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>std::list&lt; WeakSingleton&lt; Profiler &gt;::key_t &gt;</type>
      <name>m_register</name>
      <anchorfile>structmolpro_1_1profiler_1_1WeakSingleton.html</anchorfile>
      <anchor>a4783e755ed855e3d60bdb9033b22b4f3</anchor>
      <arglist></arglist>
    </member>
    <member kind="function" static="yes">
      <type>static std::shared_ptr&lt; Object &gt;</type>
      <name>single</name>
      <anchorfile>structmolpro_1_1profiler_1_1WeakSingleton.html</anchorfile>
      <anchor>ab496a3d4beee89d114bd289157702713</anchor>
      <arglist>(const std::string &amp;key, T &amp;&amp;... constructor_args)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static std::shared_ptr&lt; Object &gt;</type>
      <name>single</name>
      <anchorfile>structmolpro_1_1profiler_1_1WeakSingleton.html</anchorfile>
      <anchor>afb8cd7c6d036fdddae66cc5072ea59f3</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>erase</name>
      <anchorfile>structmolpro_1_1profiler_1_1WeakSingleton.html</anchorfile>
      <anchor>a0b23f3004d315673855ae374455d3fda</anchor>
      <arglist>(Object *obj)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>erase</name>
      <anchorfile>structmolpro_1_1profiler_1_1WeakSingleton.html</anchorfile>
      <anchor>a2c1bd92198f29ae39c2adc10bbf7c653</anchor>
      <arglist>(const std::string &amp;key)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>clear</name>
      <anchorfile>structmolpro_1_1profiler_1_1WeakSingleton.html</anchorfile>
      <anchor>af9c671c590eb3c5ab7770618447f7d23</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static std::list&lt; key_t &gt;</type>
      <name>m_register</name>
      <anchorfile>structmolpro_1_1profiler_1_1WeakSingleton.html</anchorfile>
      <anchor>a20642aa01b8a73ab3f2ee1457eecbe51</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="namespace">
    <name>molpro</name>
    <filename>namespacemolpro.html</filename>
    <namespace>molpro::profiler</namespace>
    <member kind="typedef">
      <type>profiler::Profiler</type>
      <name>Profiler</name>
      <anchorfile>namespacemolpro.html</anchorfile>
      <anchor>ac79b0d9202cceb3c8b8be0533e7ebe93</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="namespace">
    <name>molpro::profiler</name>
    <filename>namespacemolpro_1_1profiler.html</filename>
    <namespace>molpro::profiler::detail</namespace>
    <namespace>molpro::profiler::dotgraph</namespace>
    <class kind="class">molpro::profiler::Counter</class>
    <class kind="class">molpro::profiler::Node</class>
    <class kind="struct">molpro::profiler::NodePathError</class>
    <class kind="class">molpro::profiler::Profiler</class>
    <class kind="class">molpro::profiler::Timer</class>
    <class kind="struct">molpro::profiler::WeakSingleton</class>
    <member kind="enumeration">
      <type></type>
      <name>SortBy</name>
      <anchorfile>namespacemolpro_1_1profiler.html</anchorfile>
      <anchor>ae9e00894b27b784dadfc76b56d63a602</anchor>
      <arglist></arglist>
      <enumvalue file="namespacemolpro_1_1profiler.html" anchor="ae9e00894b27b784dadfc76b56d63a602a2d86bdac01a3315b95794ffa7360edc3">wall</enumvalue>
      <enumvalue file="namespacemolpro_1_1profiler.html" anchor="ae9e00894b27b784dadfc76b56d63a602ad9747e2da342bdb995f6389533ad1a3d">cpu</enumvalue>
      <enumvalue file="namespacemolpro_1_1profiler.html" anchor="ae9e00894b27b784dadfc76b56d63a602af2bb91e8e7436eaa944c378d44066a79">calls</enumvalue>
      <enumvalue file="namespacemolpro_1_1profiler.html" anchor="ae9e00894b27b784dadfc76b56d63a602aba19a09a68a66f8ad972ef8a5fba6f0d">operations</enumvalue>
      <enumvalue file="namespacemolpro_1_1profiler.html" anchor="ae9e00894b27b784dadfc76b56d63a602a334c4a4c42fdb79d7ebc3e73b517e6f8">none</enumvalue>
    </member>
    <member kind="function">
      <type>void</type>
      <name>report</name>
      <anchorfile>namespacemolpro_1_1profiler.html</anchorfile>
      <anchor>ac11f7b0c852f3f13febeb7abcf47e76c</anchor>
      <arglist>(const std::shared_ptr&lt; Node&lt; Counter &gt;&gt; &amp;root, const std::string &amp;description, std::ostream &amp;out, bool cumulative=true, SortBy sort_by=SortBy::wall)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>report</name>
      <anchorfile>namespacemolpro_1_1profiler.html</anchorfile>
      <anchor>a95116b9c056810080c206e1315c57597</anchor>
      <arglist>(const Profiler &amp;prof, std::ostream &amp;out, bool cumulative=true, SortBy sort_by=SortBy::wall)</arglist>
    </member>
    <member kind="function">
      <type>std::string</type>
      <name>get_dotgraph</name>
      <anchorfile>namespacemolpro_1_1profiler.html</anchorfile>
      <anchor>ae3e6f4ce4f648798374f80e9b3bc9de3</anchor>
      <arglist>(const Profiler &amp;prof, int hot[3], int cool[3], double threshold, bool get_percentage_time)</arglist>
    </member>
    <member kind="function">
      <type>std::ostream &amp;</type>
      <name>operator&lt;&lt;</name>
      <anchorfile>namespacemolpro_1_1profiler.html</anchorfile>
      <anchor>a14b9a5bef664e269695fc781a7243194</anchor>
      <arglist>(std::ostream &amp;os, const Profiler &amp;obj)</arglist>
    </member>
    <member kind="function">
      <type>MPI_Comm</type>
      <name>comm_self</name>
      <anchorfile>namespacemolpro_1_1profiler.html</anchorfile>
      <anchor>a33e5fcf11ef48517c70045c34eb845bc</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>MPI_Comm</type>
      <name>comm_global</name>
      <anchorfile>namespacemolpro_1_1profiler.html</anchorfile>
      <anchor>a4d67b0e0440a02ee2e44f62c44704eb6</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>size_global</name>
      <anchorfile>namespacemolpro_1_1profiler.html</anchorfile>
      <anchor>a9b51b64942c4a95951e0b462e7f95f02</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>rank_global</name>
      <anchorfile>namespacemolpro_1_1profiler.html</anchorfile>
      <anchor>a12ec6404997c63364ce0f8385e2d359c</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>init</name>
      <anchorfile>namespacemolpro_1_1profiler.html</anchorfile>
      <anchor>a65334f2187f52342168b35a4a39e276c</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>finalize</name>
      <anchorfile>namespacemolpro_1_1profiler.html</anchorfile>
      <anchor>a0f84b251e1aa8ee01afc7a3bf638ec80</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="namespace">
    <name>molpro::profiler::detail</name>
    <filename>namespacemolpro_1_1profiler_1_1detail.html</filename>
    <class kind="struct">molpro::profiler::detail::AccessCalls</class>
    <class kind="struct">molpro::profiler::detail::AccessCPU</class>
    <class kind="struct">molpro::profiler::detail::AccessOperations</class>
    <class kind="struct">molpro::profiler::detail::AccessWall</class>
    <class kind="struct">molpro::profiler::detail::Compare</class>
    <class kind="struct">molpro::profiler::detail::None</class>
    <class kind="struct">molpro::profiler::detail::TreePath</class>
    <member kind="function">
      <type>std::list&lt; std::string &gt;</type>
      <name>path_to_node</name>
      <anchorfile>namespacemolpro_1_1profiler_1_1detail.html</anchorfile>
      <anchor>ab7d6d9576c05918cef5541b9ddaf9af5</anchor>
      <arglist>(std::shared_ptr&lt; Node&lt; Counter &gt;&gt; node)</arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>total_operation_count</name>
      <anchorfile>namespacemolpro_1_1profiler_1_1detail.html</anchorfile>
      <anchor>aad8e9f89503715763086231dc2797cf3</anchor>
      <arglist>(const std::shared_ptr&lt; Node&lt; Counter &gt;&gt; &amp;node)</arglist>
    </member>
    <member kind="function">
      <type>std::string</type>
      <name>format_path_cumulative</name>
      <anchorfile>namespacemolpro_1_1profiler_1_1detail.html</anchorfile>
      <anchor>a5caed20c322a633098aebc3f36ac7d67</anchor>
      <arglist>(const std::list&lt; std::string &gt; &amp;path)</arglist>
    </member>
    <member kind="function">
      <type>std::string</type>
      <name>format_path_not_cumulative</name>
      <anchorfile>namespacemolpro_1_1profiler_1_1detail.html</anchorfile>
      <anchor>a972122028a4b2ae6f71177c7b65ebec5</anchor>
      <arglist>(const std::list&lt; std::string &gt; &amp;path)</arglist>
    </member>
    <member kind="function">
      <type>std::string</type>
      <name>format_single_path</name>
      <anchorfile>namespacemolpro_1_1profiler_1_1detail.html</anchorfile>
      <anchor>a842ead03ba3bffc6f2cbbda82b26b46c</anchor>
      <arglist>(const std::list&lt; std::string &gt; &amp;path, bool cumulative)</arglist>
    </member>
    <member kind="function">
      <type>std::map&lt; TreePath, std::shared_ptr&lt; Node&lt; Counter &gt; &gt;, CompareTreePaths &gt;</type>
      <name>sort_children</name>
      <anchorfile>namespacemolpro_1_1profiler_1_1detail.html</anchorfile>
      <anchor>a158d44cef39d6627e0ba73f34f821c25</anchor>
      <arglist>(const std::shared_ptr&lt; Node&lt; Counter &gt;&gt; &amp;root, bool cumulative)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>format_paths</name>
      <anchorfile>namespacemolpro_1_1profiler_1_1detail.html</anchorfile>
      <anchor>a88a87b883674999c31a099aff7b6adc4</anchor>
      <arglist>(std::list&lt; std::string &gt; &amp;path_names, bool append)</arglist>
    </member>
    <member kind="function">
      <type>std::string</type>
      <name>frequency</name>
      <anchorfile>namespacemolpro_1_1profiler_1_1detail.html</anchorfile>
      <anchor>a69c868ffbce95a5e750df756a4d4e710</anchor>
      <arglist>(size_t n_op, double time)</arglist>
    </member>
    <member kind="function">
      <type>std::string</type>
      <name>seconds</name>
      <anchorfile>namespacemolpro_1_1profiler_1_1detail.html</anchorfile>
      <anchor>afff5a46cc57c89e470b67badde79f97a</anchor>
      <arglist>(double time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>write_timing</name>
      <anchorfile>namespacemolpro_1_1profiler_1_1detail.html</anchorfile>
      <anchor>a6ddb1fe00da5f79b9535edb00cd0c44f</anchor>
      <arglist>(std::ostream &amp;out, double time, size_t n_op)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>write_report</name>
      <anchorfile>namespacemolpro_1_1profiler_1_1detail.html</anchorfile>
      <anchor>a4215518e44648acb4ea0d4bed1de3e5d</anchor>
      <arglist>(const Node&lt; Counter &gt; &amp;root, const std::string &amp;description, const std::list&lt; TreePath &gt; &amp;paths, std::ostream &amp;out, bool cumulative)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>write_dotgraph</name>
      <anchorfile>namespacemolpro_1_1profiler_1_1detail.html</anchorfile>
      <anchor>a181c763040422c59d6a4363c52dc50ff</anchor>
      <arglist>(std::string path, const std::string &amp;dotgraph)</arglist>
    </member>
    <member kind="function">
      <type>template std::map&lt; TreePath, std::shared_ptr&lt; Node&lt; Counter &gt; &gt;, Compare&lt; AccessWall &gt; &gt;</type>
      <name>sort_children&lt; Compare&lt; AccessWall &gt; &gt;</name>
      <anchorfile>namespacemolpro_1_1profiler_1_1detail.html</anchorfile>
      <anchor>a92f349a1cdcfae03d5fdf9b0514e6d01</anchor>
      <arglist>(const std::shared_ptr&lt; Node&lt; Counter &gt;&gt; &amp;root, bool cumulative)</arglist>
    </member>
    <member kind="function">
      <type>template std::map&lt; TreePath, std::shared_ptr&lt; Node&lt; Counter &gt; &gt;, Compare&lt; AccessCPU &gt; &gt;</type>
      <name>sort_children&lt; Compare&lt; AccessCPU &gt; &gt;</name>
      <anchorfile>namespacemolpro_1_1profiler_1_1detail.html</anchorfile>
      <anchor>a73278d1a7e33f09eab7edc0ecc009a22</anchor>
      <arglist>(const std::shared_ptr&lt; Node&lt; Counter &gt;&gt; &amp;root, bool cumulative)</arglist>
    </member>
    <member kind="function">
      <type>template std::map&lt; TreePath, std::shared_ptr&lt; Node&lt; Counter &gt; &gt;, Compare&lt; AccessCalls &gt; &gt;</type>
      <name>sort_children&lt; Compare&lt; AccessCalls &gt; &gt;</name>
      <anchorfile>namespacemolpro_1_1profiler_1_1detail.html</anchorfile>
      <anchor>a0b9fd82385051751615288e9228196e4</anchor>
      <arglist>(const std::shared_ptr&lt; Node&lt; Counter &gt;&gt; &amp;root, bool cumulative)</arglist>
    </member>
    <member kind="function">
      <type>template std::map&lt; TreePath, std::shared_ptr&lt; Node&lt; Counter &gt; &gt;, Compare&lt; AccessOperations &gt; &gt;</type>
      <name>sort_children&lt; Compare&lt; AccessOperations &gt; &gt;</name>
      <anchorfile>namespacemolpro_1_1profiler_1_1detail.html</anchorfile>
      <anchor>a9ebf3f2da57f9cae509c166d26f47870</anchor>
      <arglist>(const std::shared_ptr&lt; Node&lt; Counter &gt;&gt; &amp;root, bool cumulative)</arglist>
    </member>
    <member kind="function">
      <type>template std::map&lt; TreePath, std::shared_ptr&lt; Node&lt; Counter &gt; &gt;, Compare&lt; None &gt; &gt;</type>
      <name>sort_children&lt; Compare&lt; None &gt; &gt;</name>
      <anchorfile>namespacemolpro_1_1profiler_1_1detail.html</anchorfile>
      <anchor>a68c49ede9bd6bac29cd5d2c8607d1d1a</anchor>
      <arglist>(const std::shared_ptr&lt; Node&lt; Counter &gt;&gt; &amp;root, bool cumulative)</arglist>
    </member>
  </compound>
  <compound kind="namespace">
    <name>molpro::profiler::dotgraph</name>
    <filename>namespacemolpro_1_1profiler_1_1dotgraph.html</filename>
    <class kind="class">molpro::profiler::dotgraph::GraphEntry</class>
    <member kind="enumeration">
      <type></type>
      <name>EntryType</name>
      <anchorfile>namespacemolpro_1_1profiler_1_1dotgraph.html</anchorfile>
      <anchor>aa0b0063e1f074ee13a020772f9be9278</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>node</name>
      <anchorfile>namespacemolpro_1_1profiler_1_1dotgraph.html</anchorfile>
      <anchor>aa0b0063e1f074ee13a020772f9be9278a6359ecc7079bd086646217fa723b2953</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>edge</name>
      <anchorfile>namespacemolpro_1_1profiler_1_1dotgraph.html</anchorfile>
      <anchor>aa0b0063e1f074ee13a020772f9be9278a2890de0e5d7bbd250952da263ec702d7</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>std::string</type>
      <name>blend_colours</name>
      <anchorfile>namespacemolpro_1_1profiler_1_1dotgraph.html</anchorfile>
      <anchor>acfcdf6480f3ea96c26871d33267ef276</anchor>
      <arglist>(double ratio, int hot_colour[3], int cool_colour[3])</arglist>
    </member>
    <member kind="function">
      <type>std::string</type>
      <name>print_time</name>
      <anchorfile>namespacemolpro_1_1profiler_1_1dotgraph.html</anchorfile>
      <anchor>aa4ce68bf98b2c6aa09bcc3bbde0ebc91</anchor>
      <arglist>(double time, double total_time, bool show_percentage_time)</arglist>
    </member>
    <member kind="function">
      <type>std::string</type>
      <name>make_box</name>
      <anchorfile>namespacemolpro_1_1profiler_1_1dotgraph.html</anchorfile>
      <anchor>ab489134bfbfa7eccb77cad75a366256f</anchor>
      <arglist>(std::string name, double time, double total_time, size_t call_count, size_t opcount, int hot[3], int cool[3], bool show_percentage_time)</arglist>
    </member>
    <member kind="function">
      <type>std::string</type>
      <name>make_arrow</name>
      <anchorfile>namespacemolpro_1_1profiler_1_1dotgraph.html</anchorfile>
      <anchor>a108f2fadbabcb2f32e3dcab542e9d6c8</anchor>
      <arglist>(std::string name_from, std::string name_to, double time, double total_time, size_t call_count, int hot[3], int cool[3], bool show_percentage_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>combine_graph_entries</name>
      <anchorfile>namespacemolpro_1_1profiler_1_1dotgraph.html</anchorfile>
      <anchor>a94096778467dab3841684179dcb57e89</anchor>
      <arglist>(GraphEntry &amp;entry1, GraphEntry &amp;entry2)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>merge_vec</name>
      <anchorfile>namespacemolpro_1_1profiler_1_1dotgraph.html</anchorfile>
      <anchor>afa5aa824817491b9fb30a8ca99e1456e</anchor>
      <arglist>(std::vector&lt; GraphEntry &gt; &amp;graph_entries)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>apply_threshold</name>
      <anchorfile>namespacemolpro_1_1profiler_1_1dotgraph.html</anchorfile>
      <anchor>a7c71ad64512065278cdfbbe6956dbc86</anchor>
      <arglist>(std::vector&lt; GraphEntry &gt; &amp;graph_entries, double threshold, double total_time)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>has_parent</name>
      <anchorfile>namespacemolpro_1_1profiler_1_1dotgraph.html</anchorfile>
      <anchor>a3642a4c9eb27db3c082b371b3edfed75</anchor>
      <arglist>(GraphEntry &amp;child, std::vector&lt; GraphEntry &gt; &amp;graph_entries)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>destroy_orphans</name>
      <anchorfile>namespacemolpro_1_1profiler_1_1dotgraph.html</anchorfile>
      <anchor>abb9e15165e18618899de82ff257f6c5e</anchor>
      <arglist>(std::vector&lt; GraphEntry &gt; &amp;graph_entries)</arglist>
    </member>
    <member kind="function">
      <type>std::string</type>
      <name>get_graph_markup</name>
      <anchorfile>namespacemolpro_1_1profiler_1_1dotgraph.html</anchorfile>
      <anchor>ad8971d6e31c8eb11433fcfdb4da5f925</anchor>
      <arglist>(std::vector&lt; GraphEntry &gt; &amp;graph_entries, double total_time, int hot[3], int cool[3], bool show_percentage_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>make_dotgraph_vec</name>
      <anchorfile>namespacemolpro_1_1profiler_1_1dotgraph.html</anchorfile>
      <anchor>ad31aac833a9bffc708d5bc12e9f300f0</anchor>
      <arglist>(std::shared_ptr&lt; Node&lt; Counter &gt;&gt; root, double total_time, std::vector&lt; GraphEntry &gt; &amp;graph_entries)</arglist>
    </member>
    <member kind="function">
      <type>std::string</type>
      <name>make_dotgraph</name>
      <anchorfile>namespacemolpro_1_1profiler_1_1dotgraph.html</anchorfile>
      <anchor>a35bc62b85ffbd315947fb019ca03cb2d</anchor>
      <arglist>(std::shared_ptr&lt; Node&lt; Counter &gt;&gt; root, double total_time, int hot[3], int cool[3], double threshold, bool show_percentage_time)</arglist>
    </member>
  </compound>
  <compound kind="namespace">
    <name>profilerf</name>
    <filename>namespaceprofilerf.html</filename>
    <class kind="interface">profilerf::profiler</class>
    <member kind="function">
      <type>type(profiler) function</type>
      <name>profilernewf</name>
      <anchorfile>namespaceprofilerf.html</anchorfile>
      <anchor>a4305fc2421b258688fbf25f9990eef45</anchor>
      <arglist>(name, sort, level, comm, cpu)</arglist>
    </member>
    <member kind="function">
      <type>subroutine</type>
      <name>profilerstartf</name>
      <anchorfile>namespaceprofilerf.html</anchorfile>
      <anchor>a21bb07e910aa5210f417a96b52860dfa</anchor>
      <arglist>(this, name)</arglist>
    </member>
    <member kind="function">
      <type>subroutine</type>
      <name>profileractivef</name>
      <anchorfile>namespaceprofilerf.html</anchorfile>
      <anchor>aee6bc20402f95b8f4a39a50146e6a872</anchor>
      <arglist>(this, level, stopPrint)</arglist>
    </member>
    <member kind="function">
      <type>subroutine</type>
      <name>profilerstopf</name>
      <anchorfile>namespaceprofilerf.html</anchorfile>
      <anchor>a3c76c1bffda627c4154fef4057c9f2cc</anchor>
      <arglist>(this, name, operations)</arglist>
    </member>
    <member kind="function">
      <type>subroutine</type>
      <name>profilerprintf</name>
      <anchorfile>namespaceprofilerf.html</anchorfile>
      <anchor>a73cd46ae3b4803b7f4ead5b6566c5257</anchor>
      <arglist>(this, unit, verbosity, cumulative, precision)</arglist>
    </member>
    <member kind="function">
      <type>subroutine</type>
      <name>profilerdotgraphf</name>
      <anchorfile>namespaceprofilerf.html</anchorfile>
      <anchor>ac596773d5e35037831f6b1d881de78a8</anchor>
      <arglist>(this, path, threshold, cumulative)</arglist>
    </member>
    <member kind="function">
      <type>subroutine</type>
      <name>profilerdestroyf</name>
      <anchorfile>namespaceprofilerf.html</anchorfile>
      <anchor>a86f8942d3cef5b0770394e166a5a3176</anchor>
      <arglist>(this)</arglist>
    </member>
    <member kind="function">
      <type>subroutine</type>
      <name>destructor</name>
      <anchorfile>namespaceprofilerf.html</anchorfile>
      <anchor>a1f0cc0bdeea4ec612876877099a5a335</anchor>
      <arglist>(this)</arglist>
    </member>
    <member kind="variable">
      <type>integer, parameter, public</type>
      <name>mpicomm_kind</name>
      <anchorfile>namespaceprofilerf.html</anchorfile>
      <anchor>ae2af48b1d9953f76238633b906291c04</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="page">
    <name>index</name>
    <title>Welcome to Profiler!</title>
    <filename>index</filename>
  </compound>
</tagfile>
